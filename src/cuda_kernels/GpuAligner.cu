#include <iostream>
#include <cuda.h>
#include "GpuAligner.h"
#include <vector>
#include "nanopolish_profile_hmm_r9.h"

#define MAX_STATES 1024

__device__ float logsumexpf(float x, float y){
    if(x == -INFINITY && y == -INFINITY){
        return -INFINITY;
    }
    float result = fmax(x, y) + log1pf(expf(-fabsf(y - x)));
    return result;
}

__device__ float lp_match_r9(int rank,
                             float mean,
                             float * poreModelLevelLogStdv,
                             float * poreModelLevelStdv,
                             float * poreModelLevelMean,
                             float scale,
                             float shift,
                             float var,
                             float logVar,
                             bool debug = false){

    float log_inv_sqrt_2pi = log(0.3989422804014327); // no need to calculate this every time. better solutions available..

    // STEP 1: GET DRIFT-SCALED LEVEL:
    float level = mean;
    // TODO: Apply scaling to these 3 model values as is done in the CPP implementation
    //these can just be pulled from the model

    float gaussian_mean = scale * poreModelLevelMean[rank] + shift;
    float gaussian_stdv = poreModelLevelStdv[rank] * var;
    float gaussian_log_level_stdv = poreModelLevelLogStdv[rank] + logVar;

    // Step 3: calculate log-normal PDF
    float a = (level - gaussian_mean) / gaussian_stdv; // g is the gaussian parameters

    float emission = log_inv_sqrt_2pi - gaussian_log_level_stdv + (-0.5f * a * a); // log_inv_sqrt_2pi is defined in a comment above

    return emission; // log_inv_sqrt_2pi is defined in a comment above

}

__global__ void getScores(float * eventData,
                          float * readEventsPerBase,
                          int * numRowsPerRead,
                          int * eventStarts,
                          int * eventStrides,
                          int * kmer_ranks,
                          int * kmer_ranks_rc,
                          int * eventOffsets, // Offset to use for getting an event IDX for a specific read (read obtained by block IDX)
                          float * poreModelLevelLogStdv,
                          float * poreModelLevelStdv,
                          float * poreModelLevelMean,
                          float * scaleDev,
                          float * shiftDev,
                          float * varDev,
                          float * logVarDev,
                          float * preFlankingDev,
                          float * postFlankingDev,
                          float * returnValues) {

    // Initialise the prev probability row, which is the row of the DP table
    int n_kmers = blockDim.x; // Question: How does this deal with the case where the block is bigger than the sequence, such as if one variant is a deletion?
    int n_states = n_kmers * PSR9_NUM_STATES + 2 * PSR9_NUM_STATES; // 3 blocks per kmer and then 3 each for start and end state.

    //initialise the return value
    returnValues[blockIdx.x] = -INFINITY;

    __shared__ float prevProbabilities[MAX_STATES];

    // Initialise the previous probabilities
    for (int i = 0; i < n_states - PSR9_NUM_STATES; i++) {
        prevProbabilities[i] = -INFINITY;
    }
    for (int i = n_states - PSR9_NUM_STATES; i < n_states; i++) {
        prevProbabilities[i] = 0; // Is this correct?
    }

    //Step 1: calculate transitions. For now we are going to use external params.
    int readIdx = blockIdx.x;
    float read_events_per_base = readEventsPerBase[readIdx];
    int numRows = numRowsPerRead[readIdx]; // Number of rows in this DP table.
    int e_start = eventStarts[readIdx]; // Event start for read
    int e_stride = eventStrides[readIdx];
    bool rc = false;
    if (e_stride == -1){
        rc = true;
    }
    int e_offset = eventOffsets[readIdx]; // Within the event means etc, the offset needed for this block to get a specific event

    if(blockIdx.x==2){ // read 2 is an RC read
        printf("Block IDX is %i and stride is %i\n", blockIdx.x, e_stride);
    }

    int kmerIdx = threadIdx.x;
    uint32_t rank;

    if (rc == true) {
        rank = kmer_ranks_rc[kmerIdx];
        //printf("Using an RC rank of %i\n", rank);
    }else{
        rank = kmer_ranks[kmerIdx];
    }

    float p_stay = 1 - (1 / read_events_per_base);
    float p_skip = 0.0025;
    float p_bad = 0.001;
    float p_bad_self = p_bad;
    float p_skip_self = 0.3;

    float p_mk = p_skip; // probability of not observing an event at all
    float p_mb = p_bad; // probabilty of observing a bad event
    float p_mm_self = p_stay; // probability of observing additional events from this k-mer
    float p_mm_next = 1.0f - p_mm_self - p_mk - p_mb; // normal movement from state to state

    // transitions from event split state in previous block
    float p_bb = p_bad_self;
    float p_bk, p_bm_next, p_bm_self;
    p_bk = p_bm_next = p_bm_self = (1.0f - p_bb) / 3;

    // transitions from kmer skip state in previous block
    float p_kk = p_skip_self;
    float p_km = 1.0f - p_kk;

    // We assign some transition probabilities. I believe this is correct and they don't vary by location in the sequence
    float lp_mk = log(p_mk);
    float lp_mb = log(p_mb);
    float lp_mm_self = log(p_mm_self);
    float lp_mm_next = log(p_mm_next);
    float lp_bb = log(p_bb);
    float lp_bk = log(p_bk);
    float lp_bm_next = log(p_bm_next);
    float lp_bm_self = log(p_bm_self);
    float lp_kk = log(p_kk);
    float lp_km = log(p_km);

    float lp_sm, lp_ms;
    lp_sm = lp_ms = 0.0f;

    // Start filling out the "DP table"
    // Each thread is going to work on an individual P-HMM Block
    // WRONG - need to use threadIdx & think carefully. we have one thread per block/kmer. each block has 3 states tho.
    //int kmerIdx = blockIdx.x;
    int curBlockIdx = kmerIdx + 1; // Accounts for fact that we are not working with start block.
    int prevBlockIdx = curBlockIdx -1;
    int prevBlockOffset = PSR9_NUM_STATES * prevBlockIdx;
    int curBlockOffset = PSR9_NUM_STATES * curBlockIdx;

    // the penalty is controlled by the transition probability
    float BAD_EVENT_PENALTY = 0.0f;

    float scale = scaleDev[readIdx];
    float shift = shiftDev[readIdx];
    float var = varDev[readIdx];
    float logVar = logVarDev[readIdx];

    for(int row=1; row<numRows;row++){
        // Emission probabilities
        int event_idx = e_start + (row - 1) * e_stride;
        float event_mean = eventData[e_offset + row - 1];
        float preFlank = preFlankingDev[e_offset + row - 1];
        float postFlank = postFlankingDev[e_offset + row - 1];

        bool debug = false;

        if (threadIdx.x == 0 && (row == numRows -1) && blockIdx.x == 2){
            debug = true;
        }

        float lp_emission_m = lp_match_r9(rank,
                                          event_mean,
                                          poreModelLevelLogStdv,
                                          poreModelLevelStdv,
                                          poreModelLevelMean,
                                          scale,
                                          shift,
                                          var,
                                          logVar,
                                          debug);


        float lp_emission_b = BAD_EVENT_PENALTY;

        // Get all the scores for a match
        float HMT_FROM_SAME_M = lp_mm_self + prevProbabilities[curBlockOffset + PSR9_MATCH];
        float HMT_FROM_PREV_M = lp_mm_next + prevProbabilities[prevBlockOffset + PSR9_MATCH];
        float HMT_FROM_SAME_B = lp_bm_self + prevProbabilities[curBlockOffset + PSR9_BAD_EVENT];
        float HMT_FROM_PREV_B = lp_bm_next + prevProbabilities[prevBlockOffset + PSR9_BAD_EVENT];
        float HMT_FROM_PREV_K = lp_km + prevProbabilities[prevBlockOffset + PSR9_KMER_SKIP];



        // m_s is the probability of going from the start state
        // to this kmer. The start state is (currently) only
        // allowed to go to the first kmer. If ALLOW_PRE_CLIP
        // is defined, we allow all events before this one to be skipped,
        // with a penalty;
        float HMT_FROM_SOFT = (kmerIdx == 0 &&
                               (event_idx == e_start ||
                                (HAF_ALLOW_PRE_CLIP)))  ? lp_sm  + preFlank : -INFINITY; // TODO: Add flag for HAF ALLOW_PRE_CLIP

        if (blockIdx.x == 2 && threadIdx.x == 0 && row == 2){
            printf("HMT_FROM_SOFT should be (?) -5.99 but is in fact %f\n", HMT_FROM_SOFT);
            printf("event IDX is %i\n", event_idx);
            printf("e_start is %i\n", e_start);
        }

        // calculate the score
        float sum = HMT_FROM_SAME_M;
        sum = logsumexpf(sum, HMT_FROM_SOFT);
        sum = logsumexpf(sum, HMT_FROM_PREV_M);
        sum = logsumexpf(sum, HMT_FROM_SAME_B);
        sum = logsumexpf(sum, HMT_FROM_PREV_B);
        sum = logsumexpf(sum, HMT_FROM_PREV_K);
        sum += lp_emission_m;


        float newMatchScore = sum;
        // Here need to calculate the bad event score


        // state PSR9_BAD_EVENT
        HMT_FROM_SAME_M = lp_mb + prevProbabilities[curBlockOffset + PSR9_MATCH];
        HMT_FROM_PREV_M = -INFINITY; // not allowed
        HMT_FROM_SAME_B = lp_bb + prevProbabilities[curBlockOffset + PSR9_BAD_EVENT];
        HMT_FROM_PREV_B = -INFINITY;
        HMT_FROM_PREV_K = -INFINITY;
        HMT_FROM_SOFT = -INFINITY;

        sum = HMT_FROM_SAME_M;
        sum = logsumexpf(sum, HMT_FROM_PREV_M);
        sum = logsumexpf(sum, HMT_FROM_SAME_B);
        sum = logsumexpf(sum, HMT_FROM_PREV_B);
        sum = logsumexpf(sum, HMT_FROM_PREV_K);
        sum = logsumexpf(sum, HMT_FROM_SOFT);
        sum += lp_emission_b;

        float newBadEventScore = sum;

        // Write row out. prevProbabilities now becomes "current probabilities" for evaluating skips.
        prevProbabilities[curBlockOffset + PSR9_MATCH] = newMatchScore;
        prevProbabilities[curBlockOffset + PSR9_BAD_EVENT] = newBadEventScore;
        __syncthreads();

        // state PSR9_KMER_SKIP
        HMT_FROM_SAME_M = -INFINITY;
        HMT_FROM_PREV_M = lp_mk + prevProbabilities[prevBlockOffset + PSR9_MATCH];
        HMT_FROM_SAME_B = -INFINITY;
        HMT_FROM_PREV_B = lp_bk + prevProbabilities[prevBlockOffset + PSR9_BAD_EVENT];
        HMT_FROM_SOFT = -INFINITY;

        sum = HMT_FROM_SAME_M;
        sum = logsumexpf(sum, HMT_FROM_PREV_M);
        sum = logsumexpf(sum, HMT_FROM_SAME_B);
        sum = logsumexpf(sum, HMT_FROM_PREV_B);
        sum = logsumexpf(sum, HMT_FROM_PREV_K);
        sum = logsumexpf(sum, HMT_FROM_SOFT);
        sum += 0.0;//No emission. redundant.

        float newSkipScore = sum;

        prevProbabilities[curBlockOffset + PSR9_KMER_SKIP] = newSkipScore;
        __syncthreads();

        //Now need to do the skip-skip transition, which is serial so for now letting one thread execute it.
        if (threadIdx.x == 0){
            for (int blkidx = 2;blkidx <= blockDim.x; blkidx++){
                auto skipIdx = blkidx * PSR9_NUM_STATES + PSR9_KMER_SKIP;
                float prevSkipScore = prevProbabilities[skipIdx - PSR9_NUM_STATES];
                float curSkipScore = prevProbabilities[skipIdx];
                HMT_FROM_PREV_K = lp_kk + prevSkipScore;
                newSkipScore = logsumexpf(curSkipScore, HMT_FROM_PREV_K);
                prevProbabilities[skipIdx] = newSkipScore;
                __syncthreads();
            }
        }

        __syncthreads();

        int lastKmerIdx = n_kmers -1;
        int lastRowIdx = numRows -1;
        float end;
        // Now do the post-clip transition
        if(kmerIdx == lastKmerIdx && ( (HAF_ALLOW_POST_CLIP) || row == lastRowIdx)) {
            float lp1 = lp_ms + prevProbabilities[curBlockOffset + PSR9_MATCH] + postFlank;
            float lp2 = lp_ms + prevProbabilities[curBlockOffset + PSR9_BAD_EVENT] + postFlank;
            float lp3 = lp_ms + prevProbabilities[curBlockOffset + PSR9_KMER_SKIP] + postFlank;
//
//            printf(">GPU Post-clip transition on row %i, read %i, threadIdx is %i\n"
//                           "LP1=%f\n"
//                           "LP2=%f\n"
//                           "LP3=%f\n",
//                   row,
//                   blockIdx.x,
//                   threadIdx.x,
//                   lp1,
//                   lp2,
//                   lp3);

            end = returnValues[blockIdx.x];
            end = logsumexpf(end, lp1);
            end = logsumexpf(end, lp2);
            end = logsumexpf(end, lp3);
            returnValues[blockIdx.x] = end;
        }
        // Now do the end state
        __syncthreads();

        if ((blockIdx.x == 2) && (threadIdx.x == 0)){
//            printf("rank %i\n", rank);
//            printf("event mean %f\n", event_mean);
//            printf("poreModelLevelLogStdv %f\n", poreModelLevelLogStdv[0]);
//            printf("poreModelLevelStdv %f\n", poreModelLevelStdv[0]);
//            printf("poreModelLevelMean %f\n", poreModelLevelMean[0]);
//            printf("lp_emission_m is %f\n", lp_emission_m);
//            printf("PSR9_MATCH is %i\n", PSR9_MATCH);
//            printf(">GPU score HMT_FROM_SAME_M is %f\n", HMT_FROM_SAME_M);
//            printf(">GPU score HMT_FROM_PREV_M is %f\n", HMT_FROM_PREV_M);
//            printf(">GPU score HMT_FROM_SAME_B is %f\n", HMT_FROM_SAME_B);
//            printf(">GPU score HMT_FROM_PREV_B is %f\n", HMT_FROM_PREV_B);
//            printf(">GPU score HMT_FROM_PREV_K is %f\n", HMT_FROM_PREV_K);
//            printf(">GPU newSkipScore is %f\n", newSkipScore);
//            printf("Number of states is %i\n", n_states);
                for (int c = 0; c < n_states; c++) {
                    printf("GPU> Value for row %i and col %i is %f\n", row, c, prevProbabilities[c]);
                }
            printf("HMT_FROM_SOFT = %f\n", HMT_FROM_SOFT);
            printf("lp_mk = %f\n", lp_mk);
            printf("lp_mb = %f\n", lp_mb);
            printf("lp_mm_self = %f\n", lp_mm_self);
            printf("lp_mm_next = %f\n", lp_mm_next);
            printf("lp_bb = %f\n", lp_bb);
            printf("lp_bk = %f\n", lp_bk);
            printf("lp_bm_next = %f\n", lp_bm_next);
            printf("lp_bm_self = %f\n", lp_bm_self);
            printf("lp_kk = %f\n", lp_kk);
            printf("lp_km = %f\n", lp_km);

        }
        }
        __syncthreads();
}


GpuAligner::GpuAligner()
{
    y = 20;
    asize = y*sizeof(int);
    for (int i=0; i<y; i++)
        n[i] = i;
}

std::vector<double> scoreKernel(std::vector<HMMInputSequence> sequences,
                    std::vector<HMMInputData> event_sequences,
                    uint32_t alignment_flags){

    // Extract the pore model.
    //Let's assume that every event sequence has the same pore model
    //event_sequences[0].pore_model.

    int num_reads = event_sequences.size();
    // These asserts are here during the development phase
    assert(!sequences.empty());
    assert(std::string(sequences[0].get_alphabet()->get_name()) == "nucleotide");
    for (auto e: event_sequences) {
        assert(std::string(e.pore_model->pmalphabet->get_name()) == "nucleotide");
        assert(e.read->pore_type == PT_R9);
        assert( (e.rc && e.event_stride == -1) || (!e.rc && e.event_stride == 1));
    }

    size_t num_models = sequences.size();
    double num_model_penalty = log(num_models);

    assert(num_models == 1); //this is temporary

    auto sequence = sequences[0]; // temporary. We are only going to score one sequence against a set of events for now.

    const uint32_t k = event_sequences[0].pore_model->k; //k is the kmerity
    uint32_t n_kmers = sequence.length() - k + 1; //number of kmers in the sequence

    uint32_t n_states = PSR9_NUM_STATES * (n_kmers + 2); // + 2 for explicit terminal states

    std::vector<uint32_t> n_rows; //number of rows in the DP table (n_events + 1)
    std::vector<uint32_t> e_starts; //event starts
    std::vector<int> event_strides;

    std::vector<std::vector<float>> pre_flanks;
    std::vector<std::vector<float>> post_flanks;

    for(auto e: event_sequences){
        uint32_t e_start = e.event_start_idx;
        e_starts.push_back(e_start);

        uint32_t e_stride = e.event_stride;
        event_strides.push_back(e_stride);

        uint32_t e_end = e.event_stop_idx;
        uint32_t n_events = 0;
        if(e_end > e_start)
            n_events = e_end - e_start + 1;
        else
            n_events = e_start - e_end + 1;

        n_rows.push_back(n_events + 1);

        std::vector<float> pre_flank = make_pre_flanking(e, e_start, n_events);
        std::vector<float> post_flank = make_post_flanking(e, e_start, n_events);

        pre_flanks.push_back(pre_flank);
        post_flanks.push_back(post_flank);
    }

    std::vector<uint32_t> kmer_ranks(n_kmers);
    std::vector<uint32_t> kmer_ranks_rc(n_kmers);
    for(size_t ki = 0; ki < n_kmers; ++ki) {
        kmer_ranks[ki] = sequences[0].get_kmer_rank(ki, k, false);
        kmer_ranks_rc[ki] = sequences[0].get_kmer_rank(ki, k, true);
    }

    // Prepare raw data and send it over to the score calculator kernel

    // Buffer 1: Raw event data and associated starts and stops

    size_t numEventsTotal = 0;
    //1. Count the total number of events across all reads
    std::vector<int> eventLengths;
    std::vector<float> eventsPerBase;
    for (auto e: event_sequences){
        size_t numEvents = e.read->events->size();
        float readEventsPerBase = e.read->events_per_base[e.strand];

        eventLengths.push_back(numEvents);
        eventsPerBase.push_back(readEventsPerBase);

        numEventsTotal += numEvents;
    }


    //Allocate a host buffer to store the event means, pre and post-flank data
    float * eventMeans;
    size_t eventMeansSize = numEventsTotal * sizeof(float);
    cudaHostAlloc(&eventMeans, eventMeansSize , cudaHostAllocDefault);

    //Allocate a host buffer to store the event means, pre and post-flank data
    float * preFlankingHost;
    float * postFlankingHost;
    cudaHostAlloc(&preFlankingHost, numEventsTotal * sizeof(float) , cudaHostAllocDefault);
    cudaHostAlloc(&postFlankingHost, numEventsTotal * sizeof(float) , cudaHostAllocDefault);

    std::vector<int> eventOffsets;
    size_t offset = 0;
    for(int j=0;j<event_sequences.size();j++){
        auto ev = event_sequences[j];
        eventOffsets.push_back(offset);
        size_t num_events = 100;//TODO: FIX! ev.read->events->size();
        for (int i=0;i<num_events;i++) {
            auto event_idx =  e_starts[j] + i * event_strides[j];
            auto scaled = ev.read->get_drift_scaled_level(event_idx, ev.strand); // send the data in drift scaled
            //auto unscaled = ev.read->events[0][i].mean; //taking the first element. Not sure what the second one is..
            eventMeans[offset + i] = scaled;
            preFlankingHost[offset + i] = pre_flanks[j][i]; //also copy over the pre-flanking data, since it has a 1-1 correspondence with events
            postFlankingHost[offset + i] = post_flanks[j][i]; //also copy over the pre-flanking data, since it has a 1-1 correspondence with events
        }
        offset += num_events;
    }

    int num_states = event_sequences[0].pore_model->states.size();

    std::vector<float> pore_model_level_log_stdv(num_states);
    std::vector<float> pore_model_level_mean(num_states);
    std::vector<float> pore_model_level_stdv(num_states);

    //TODO: Fix this.
    for(int st=0; st<num_states; st++){
        auto params = event_sequences[0].pore_model->states[st]; //TODO: Is this OK?
        pore_model_level_log_stdv[st] = params.level_log_stdv;
        pore_model_level_mean[st] = params.level_mean;
        pore_model_level_stdv[st] = params.level_stdv;
    }

    std::vector<float> scale(num_reads);
    std::vector<float> shift(num_reads);
    std::vector<float> var(num_reads);
    std::vector<float> log_var(num_reads);

    for (int i=0;i<num_reads;i++){
        auto read = event_sequences[i];
        scale[i] = event_sequences[i].read->scalings[read.strand].scale;
        shift[i] = event_sequences[i].read->scalings[read.strand].shift;
        var[i] = event_sequences[i].read->scalings[read.strand].var;
        log_var[i] = event_sequences[i].read->scalings[read.strand].log_var;
    }

    float* scaleDev;
    float* shiftDev;
    float* varDev;
    float* logVarDev;

    cudaMalloc( (void**)&scaleDev, scale.size() * sizeof(float));
    cudaMalloc( (void**)&shiftDev, shift.size() * sizeof(float));
    cudaMalloc( (void**)&varDev, var.size() * sizeof(float));
    cudaMalloc( (void**)&logVarDev, log_var.size() * sizeof(float));

    cudaMemcpyAsync( scaleDev, scale.data(), scale.size() * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpyAsync( shiftDev, shift.data(), shift.size() * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpyAsync( varDev, var.data(), var.size() * sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpyAsync( logVarDev, log_var.data(), log_var.size() * sizeof(float), cudaMemcpyHostToDevice );


    float* poreModelLevelLogStdvDev;
    cudaMalloc( (void**)&poreModelLevelLogStdvDev, pore_model_level_log_stdv.size() * sizeof(float));
    cudaMemcpyAsync( poreModelLevelLogStdvDev, pore_model_level_log_stdv.data(), pore_model_level_log_stdv.size() * sizeof(float), cudaMemcpyHostToDevice );

    float* poreModelLevelMeanDev;
    cudaMalloc( (void**)&poreModelLevelMeanDev, pore_model_level_mean.size() * sizeof(float));
    cudaMemcpyAsync( poreModelLevelMeanDev, pore_model_level_mean.data(), pore_model_level_mean.size() * sizeof(float), cudaMemcpyHostToDevice );

    float* poreModelLevelStdvDev;
    cudaMalloc( (void**)&poreModelLevelStdvDev, pore_model_level_stdv.size() * sizeof(float));
    cudaMemcpyAsync( poreModelLevelStdvDev, pore_model_level_stdv.data(), pore_model_level_stdv.size() * sizeof(float), cudaMemcpyHostToDevice );


    float* eventsPerBaseDev;
    cudaMalloc( (void**)&eventsPerBaseDev, eventsPerBase.size() * sizeof(float));
    cudaMemcpyAsync( eventsPerBaseDev, eventsPerBase.data(), eventsPerBase.size() * sizeof(float), cudaMemcpyHostToDevice );

    float* eventMeansDev;
    cudaMalloc( (void**)&eventMeansDev, eventMeansSize);
    cudaMemcpyAsync( eventMeansDev, eventMeans, eventMeansSize, cudaMemcpyHostToDevice ); //malloc is taking 300us

    float* preFlankingDev;
    cudaMalloc( (void**)&preFlankingDev, eventMeansSize);
    cudaMemcpyAsync( preFlankingDev, preFlankingHost, eventMeansSize, cudaMemcpyHostToDevice ); //malloc is taking 300us

    float* postFlankingDev;
    cudaMalloc( (void**)&postFlankingDev, eventMeansSize);
    cudaMemcpyAsync( postFlankingDev, postFlankingHost, eventMeansSize, cudaMemcpyHostToDevice ); //malloc is taking 300us

    int* numRowsDev;
    cudaMalloc( (void**)&numRowsDev, n_rows.size() * sizeof(int));
    cudaMemcpyAsync( numRowsDev, n_rows.data(), n_rows.size() * sizeof(int), cudaMemcpyHostToDevice );

    int* kmerRanksDev;
    int* kmerRanksRCDev;
    cudaMalloc( (void**)&kmerRanksDev, kmer_ranks.size() * sizeof(int));
    cudaMalloc( (void**)&kmerRanksRCDev, kmer_ranks_rc.size() * sizeof(int));
    cudaMemcpyAsync( kmerRanksDev, kmer_ranks.data(), kmer_ranks.size() * sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpyAsync( kmerRanksRCDev, kmer_ranks_rc.data(), kmer_ranks_rc.size() * sizeof(int), cudaMemcpyHostToDevice );

    int* eventStartsDev;
    cudaMalloc( (void**)&eventStartsDev, e_starts.size() * sizeof(int));
    cudaMemcpyAsync( eventStartsDev, e_starts.data(), e_starts.size() * sizeof(int), cudaMemcpyHostToDevice );

    int* eventStridesDev;
    cudaMalloc( (void**)&eventStridesDev, event_strides.size() * sizeof(int));
    cudaMemcpyAsync( eventStridesDev, event_strides.data(), event_strides.size() * sizeof(int), cudaMemcpyHostToDevice );

    int* eventOffsetsDev;
    cudaMalloc( (void**)&eventOffsetsDev, eventOffsets.size() * sizeof(int));
    cudaMemcpyAsync( eventOffsetsDev, eventOffsets.data(), eventOffsets.size() * sizeof(int), cudaMemcpyHostToDevice );

    int num_blocks = n_states / PSR9_NUM_STATES;
    uint32_t num_kmers = num_blocks - 2; // two terminal blocks. Not currently used but left here for now.

    dim3 dimBlock(num_blocks - 2); // One thread per state, not including Start and Terminal state.
    dim3 dimGrid(num_reads); // let's look at only the first read

    float * returnValuesDev;
    cudaMalloc((void **) &returnValuesDev, sizeof(float) * num_reads); //one score per read

    float* returnedValues;
    cudaHostAlloc(&returnedValues, num_reads * sizeof(float) , cudaHostAllocDefault);

    printf("About to run getscores...\n");
    getScores<<<dimGrid, dimBlock>>>(eventMeansDev,
            eventsPerBaseDev,
            numRowsDev,
            eventStartsDev,
            eventStridesDev,
            kmerRanksDev,
            kmerRanksRCDev,
            eventOffsetsDev,
            poreModelLevelLogStdvDev,
            poreModelLevelStdvDev,
            poreModelLevelMeanDev,
            scaleDev,
            shiftDev,
            varDev,
            logVarDev,
            preFlankingDev,
            postFlankingDev,
            returnValuesDev);

    cudaDeviceSynchronize();
    cudaMemcpyAsync(returnedValues, returnValuesDev, num_reads *sizeof(float), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(eventMeansDev);
    cudaFree(eventsPerBaseDev);
    cudaFree(numRowsDev);
    cudaFree(eventStartsDev);
    cudaFree(eventStridesDev);
    cudaFree(kmerRanksDev);
    cudaFree(kmerRanksRCDev);
    cudaFree(eventOffsetsDev);
    cudaFree(poreModelLevelLogStdvDev);
    cudaFree(poreModelLevelStdvDev);
    cudaFree(poreModelLevelMeanDev);
    cudaFree(scaleDev);
    cudaFree(shiftDev);
    cudaFree(varDev);
    cudaFree(logVarDev);
    cudaFree(preFlankingDev);
    cudaFree(postFlankingDev);

    //Free host memory
    cudaFreeHost(eventMeans);

    //Send all the scores back
    std::vector<double> r(num_reads);
    for(int i=0; i<num_reads;i++){
        r[i]= (double) returnedValues[i];
    }

    return r;
}

std::vector<double> GpuAligner::variantScoresThresholded(std::vector<Variant> input_variants,
                                                        Haplotype base_haplotype,
                                                        std::vector<HMMInputData> event_sequences,
                                                        uint32_t alignment_flags,
                                                        int screen_score_threshold,
                                                        std::vector<std::string> methylation_types) {
    int numVariants = input_variants.size();

    std::vector<Variant> out_variants = input_variants;
    std::vector<Haplotype> variant_haplotypes(numVariants, base_haplotype);

    //loop over the vector, applying the variants to the haplotypes
    for (int i = 0; i<input_variants.size();i++){
        variant_haplotypes[i].apply_variant(input_variants[i]);
    }

    // Make methylated versions of each input sequence. Once for the base haplotype and once each for each variant
    std::vector<HMMInputSequence> base_sequences = generate_methylated_alternatives(base_haplotype.get_sequence(),
                                                                                    methylation_types);
    std::vector<std::vector<HMMInputSequence>> variant_sequences;

    for (auto v: variant_haplotypes){
        auto variant_sequence = generate_methylated_alternatives(v.get_sequence(), methylation_types);
        variant_sequences.push_back(variant_sequence);
    }

    assert(base_sequences.size() == 1);

    // return the sum of the score for the base sequences over all the event sequences
    auto base_scores = scoreKernel(base_sequences, event_sequences, alignment_flags);

    std::vector<double> v(variant_sequences.size());
    for (int i=0; i<variant_sequences.size(); i++){
        auto scores = scoreKernel(variant_sequences[i], event_sequences, alignment_flags); //TODO: Base sequence needs to be replaced with the variant itself

        double totalScore = 0.0;
        for(int k=0; k<scores.size(); k++){
            if (fabs(totalScore) < screen_score_threshold){ //threshold, hardcoded. TODO: Make this an argument, although the thresholding doesn't make much sense on GPU
                totalScore += (scores[k] - base_scores[k]);
            }
        }
        v[i] = totalScore;
    }

    return v;
}
