#include <iostream>
#include <cuda.h>
#include "GpuAligner.h"
#include <vector>
#include "nanopolish_profile_hmm_r9.h"


__global__ void findSumToN(int *n, int limit)
{
    //printf("HELLO FROM SUM\n");
    int tId = threadIdx.x;

    for (int i=0; i<=(int)log2((double)limit); i++)
    {
        if (tId%(int)(pow(2.0,(double)(i+1))) == 0){
            if (tId+(int)pow(2.0, (double)i) >= limit) break;
            n[tId] += n[tId+(int)pow(2.0, (double)i)];
        }
        __syncthreads();
    }
}


__global__ void getScores(float * eventData, float * returnValues)
{
    int tId = threadIdx.x;
    if (tId == 0) {
        printf("data: %f\n", eventData[0]);
        printf("data: %f\n", eventData[1]);
        printf("data: %f\n", eventData[2]);
    }
    returnValues[0] = 0.356;
    //__syncthreads();
}


GpuAligner::GpuAligner()
{
    y = 20;
    asize = y*sizeof(int);
    for (int i=0; i<y; i++)
        n[i] = i;
}

int GpuAligner::calculateSum()
{
    int *n_d;
    cudaMalloc( (void**)&n_d, asize );

    cudaMemcpy(n_d, n, asize, cudaMemcpyHostToDevice );

    dim3 dimBlock( y, 1 );
    dim3 dimGrid( 1, 1 );
    findSumToN<<<dimGrid, dimBlock>>>(n_d, y);
    cudaMemcpy(n, n_d, asize, cudaMemcpyDeviceToHost);
    cudaFree (n_d);
    return n[0];
}

void GpuAligner::setY(int newVal)
{
    y = newVal;
    asize = y*sizeof(int);
    for (int i=0; i<y; i++)
        n[i] = i;

}

double scoreKernel(std::vector<HMMInputSequence> sequences,
                   std::vector<HMMInputData> event_sequences,
                   uint32_t alignment_flags){

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

    for(auto e: event_sequences){
        uint32_t e_start = e.event_start_idx;
        e_starts.push_back(e_start);
        uint32_t e_end = e.event_stop_idx;
        uint32_t n_events = 0;
        if(e_end > e_start)
            n_events = e_end - e_start + 1;
        else
            n_events = e_start - e_end + 1;

        n_rows.push_back(n_events + 1);
    }


    // Prepare raw data and send it over to the score calculator kernel

    // Buffer 1: Raw event data and associated starts and stops

    size_t numEventsTotal = 0;
    //1. Count the total number of events across all reads
    std::vector<int> eventLengths;
    for (auto e: event_sequences){
        size_t numEvents = e.read->events->size();

        eventLengths.push_back(numEvents);
        numEventsTotal += numEvents;
    }

    float * eventMeans;
    //Allocate a host buffer to store the event means
    size_t eventMeansSize = numEventsTotal * sizeof(float);
    cudaHostAlloc(&eventMeans, eventMeansSize , cudaHostAllocDefault);

    size_t offset = 0;
    for (auto ev: event_sequences){
        size_t num_events = ev.read->events->size();
        for (int i=0;i<num_events;i++) {
            eventMeans[offset + i] = ev.read->events[0][i].mean; //taking the first element. Not sure what the second one is..
        }
        offset += num_events;
    }


    float* devicePtr;
    cudaMalloc( (void**)&devicePtr, eventMeansSize);
    cudaMemcpy( devicePtr, eventMeans, eventMeansSize, cudaMemcpyHostToDevice );

    dim3 dimBlock( 1, 1 );
    dim3 dimGrid( 1, 1 );

    float * returnValues;
    cudaMalloc((void **) &returnValues, sizeof(float) * num_models); //one score per read

    float * returnedValues;
    getScores<<<dimGrid, dimBlock>>>(devicePtr, returnValues);

    cudaMemcpy(returnedValues, returnValues, num_models *sizeof(float), cudaMemcpyDeviceToHost );

    auto r = returnedValues[0];
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


    //variant_haplotype.apply_variant(input_variant);

    // Make methylated versions of each input sequence. Once for the base haplotype and once each for each variant
    std::vector<HMMInputSequence> base_sequences = generate_methylated_alternatives(base_haplotype.get_sequence(),
                                                                                    methylation_types);

    assert(base_sequences.size() == 1);

    std::vector<std::vector<HMMInputSequence>> methylatedVariantSequences;
    for(auto variant: variant_haplotypes) {
        std::vector<HMMInputSequence> variant_sequences = generate_methylated_alternatives(
                variant.get_sequence(), methylation_types);
        methylatedVariantSequences.push_back(variant_sequences);

    }

    //For now let's not worry about methylation
    assert(methylatedVariantSequences.size() == numVariants);
    for (auto m: methylatedVariantSequences) {
        assert(m.size() == 1);
    }
    //Next we need to get the scores.

    // return the sum of the score for the base sequences over all the event sequences
    double base_score = scoreKernel(base_sequences, event_sequences, alignment_flags);

    std::vector<double> v;
    v.push_back(base_score);
    return v;
}