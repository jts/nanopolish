#include <iostream>
#include <cuda.h>
#include "GpuAligner.h"
#include <vector>
#include "nanopolish_profile_hmm_r9.h"

#define MAX_STATES 1024

__device__ float logsumexpf(float x, float y){
    return fmax(x, y) + log1pf(expf(-fabsf(y-x)));
}

//TODO: Implement, inc pore model
__device__ float lp_match_r9(int rank,
                             float mean,
                             float * poreModelLevelLogStdv,
                             float * poreModelLevelStdv,
                             float * poreModelLevelMean){

    float log_inv_sqrt_2pi = log(0.3989422804014327); // no need to calculate this every time. better solutions available..

    // STEP 1: GET DRIFT-SCALED LEVEL:
    float level = mean; //TODO: Do actual drift scaling. this is a cheat
    // TODO: Apply scaling to these 3 model values as is done in the CPP implementation
    //these can just be pulled from the model
    float gaussian_mean = poreModelLevelMean[rank];
    float gaussian_stdv = poreModelLevelStdv[rank];
    float gaussian_log_level_stdv = poreModelLevelLogStdv[rank];
    // Step 3: calculate log-normal PDF
    float a = (level - gaussian_mean) / gaussian_stdv; // g is the gaussian parameters
    return log_inv_sqrt_2pi - gaussian_log_level_stdv + (-0.5f * a * a); // log_inv_sqrt_2pi is defined in a comment above
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
                          float * returnValues)
{
    // kmer probabilities will be stored here
    __shared__ float prevProbabilities[MAX_STATES];
    for (int i =0;i<MAX_STATES;i++){
        prevProbabilities[i] = -INFINITY;
    }

    //float log_inv_sqrt_2pi = log(0.3989422804014327);

    //Step 1: calculate transitions. For now we are going to use external params.
    int readIdx = blockIdx.x;
    float read_events_per_base = readEventsPerBase[readIdx];
    int numRows = numRowsPerRead[readIdx]; // Number of rows in this DP table.
    int e_start = eventStarts[readIdx]; // Event start for read
    int e_stride = eventStrides[readIdx];
    int e_offset = eventOffsets[readIdx]; // Within the event means etc, the offset needed for this block to get a specific event
    //int kmer_ranks = kmerRanks[readIdx.x]; // TODO: Use RC for RC reads

    int kmerIdx = threadIdx.x;

    float p_stay = 1 - (1 / read_events_per_base);

    //printf("Events per base: %f \n", read_events_per_base);
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

    // We assign some transition probabilities. I believe this is correct and they don't vary by location in the sequence (why would they)
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

    for(int row=1; row<numRows;row++){
        // Emission probabilities
        int event_idx = e_start + (row - 1) * e_stride;
        uint32_t rank = kmer_ranks[kmerIdx]; // lexical rank of a kmer
        float event_mean = eventData[e_offset + row];
        float lp_emission_m = lp_match_r9(rank,
                                          event_mean,
                                          poreModelLevelLogStdv,
                                          poreModelLevelStdv,
                                          poreModelLevelMean);
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
        // TODO: Implemnet the HMT_FROM_SOFT score. this appears needed but I don't yet understand it.

        // NOW calculate the score
        float sum = HMT_FROM_SAME_M;
        sum = logsumexpf(sum, HMT_FROM_PREV_M);
        sum = logsumexpf(sum, HMT_FROM_SAME_B);
        sum = logsumexpf(sum, HMT_FROM_PREV_B);
        sum = logsumexpf(sum, HMT_FROM_PREV_K);
        sum += lp_emission_m;

        __syncthreads();
        prevProbabilities[curBlockIdx + PSR9_MATCH] = sum;
        __syncthreads();
    }


    returnValues[blockIdx.x] = 0.356;
    __syncthreads();
}


GpuAligner::GpuAligner()
{
    y = 20;
    asize = y*sizeof(int);
    for (int i=0; i<y; i++)
        n[i] = i;
}

double scoreKernel(std::vector<HMMInputSequence> sequences,
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
    std::vector<uint32_t> event_strides;

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

    //Allocate a host buffer to store the event means
    float * eventMeans;
    size_t eventMeansSize = numEventsTotal * sizeof(float);
    cudaHostAlloc(&eventMeans, eventMeansSize , cudaHostAllocDefault);

    std::vector<int> eventOffsets;
    size_t offset = 0;
    for (auto ev: event_sequences){
        eventOffsets.push_back(offset);
        size_t num_events = ev.read->events->size();
        for (int i=0;i<num_events;i++) {
            eventMeans[offset + i] = ev.read->events[0][i].mean; //taking the first element. Not sure what the second one is..
        }
        offset += num_events;
    }

    int num_states = event_sequences[0].pore_model->states.size();

    std::vector<float> pore_model_level_log_stdv(num_states);
    std::vector<float> pore_model_level_mean(num_states);
    std::vector<float> pore_model_level_stdv(num_states);

    for(int st=0; st<num_states; st++){
        auto params = event_sequences[0].pore_model->states[0]; //let's just initially get the params for AAAAAA
        pore_model_level_log_stdv[st] = params.level_log_stdv;
        pore_model_level_mean[st] = params.level_mean;
        pore_model_level_stdv[st] = params.level_stdv;
    }


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
    uint32_t num_kmers = num_blocks - 2; // two terminal blocks

    dim3 dimBlock(num_blocks - 2);
    dim3 dimGrid(1); // One thread per state, not including Start and Terminal state.

    float * returnValues;
    cudaMalloc((void **) &returnValues, sizeof(float) * num_reads); //one score per read

    float* returnedValues;// = new float[num_reads];
    //size_t eventMeansSize = numEventsTotal * sizeof(float);
    cudaHostAlloc(&returnedValues, num_reads * sizeof(float) , cudaHostAllocDefault);

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
            returnValues);

    //cudaDeviceSynchronize();
    cudaMemcpyAsync(returnedValues, returnValues, num_reads *sizeof(float), cudaMemcpyDeviceToHost);

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

    //Free host memory
    cudaFreeHost(eventMeans);

    float r = 0.0;
    for(int i=0; i<num_reads;i++){
        r += returnedValues[i];
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
    double base_score = scoreKernel(base_sequences, event_sequences, alignment_flags);

    std::vector<double> v(variant_sequences.size());
    for (int i=0; i<variant_sequences.size(); i++){
        double score = scoreKernel(variant_sequences[i], event_sequences, alignment_flags); //TODO: Base sequence needs to be replaced with the variant itself
        v[i] = (score - base_score);
    }

    return v;
}
