#include <iostream>
#include <cuda.h>
#include "GpuAligner.h"
#include <vector>

__global__ void findSumToN(int *n, int limit)
{
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

    assert(!sequences.empty());
    assert(std::string(sequences[0].get_alphabet()->get_name()) == "nucleotide");
    for (auto e: event_sequences) {
        assert(std::string(e.pore_model->pmalphabet->get_name()) == "nucleotide");
        assert(e.read->pore_type == PT_R9);
    }

    size_t num_models = sequences.size();
    double num_model_penalty = log(num_models);

    assert(num_models == 1); //this is temporary

    // start preparing the data for the CUDA Kernel



    return 0.210964;
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