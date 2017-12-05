//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_duration_model -- Model the duration
// of bases passing through the pore
//
#include "nanopolish_duration_model.h"
#include "nanopolish_profile_hmm.h"

DurationModel::DurationModel()
{

}

std::vector<double> DurationModel::generate_aligned_durations(const std::string& sequence,
                                                              const HMMInputData& data,
                                                              const uint32_t alignment_flags)
{
    size_t k = data.pore_model->k;
    size_t num_kmers = sequence.size() - k + 1;
    // initialize the vector of durations
    std::vector<double> duration_by_kmer_position(num_kmers, 0.0);
    std::vector<HMMAlignmentState> alignment = profile_hmm_align(sequence, data, alignment_flags);
    for(size_t ai = 0; ai < alignment.size(); ai++) {

        /*   
        fprintf(stderr, "alignment[%zu]: %s %zu %c %.5lf\n", ai, 
                                                             sequence.substr(alignment[ai].kmer_idx, 6).c_str(), 
                                                             alignment[ai].event_idx, 
                                                             alignment[ai].state, 
                                                             data.read->get_duration(alignment[ai].event_idx, data.strand));
        */
        if(alignment[ai].state != 'K') {
            duration_by_kmer_position[alignment[ai].kmer_idx] += data.read->get_duration(alignment[ai].event_idx, data.strand);
        }
    }
    return duration_by_kmer_position;
}

double DurationModel::log_gamma_sum(double x, double n)
{
    const static double a = 2.461964; // shape
    const static double b = 587.2858; // rate
    GammaParameters params;
    params.shape = a;
    params.rate = b;
    return log_gamma_sum(x, params, n);
}

double DurationModel::log_gamma_sum(double x, const GammaParameters& params, double n)
{
    assert(x >= 0.0);

    double na = n * params.shape;
    return (na * log(params.rate)) - lgamma(na) + (na - 1) * log(x) - params.rate * x;
}

GammaParameters DurationModel::gamma_fit(const std::vector<double>& input)
{
    double s = gamma_fit_calculate_s(input);
    //double k = (3 - s + sqrt(pow(s - 3.0, 2.0) + 24.0 * s)) / (12 * s);
    double k = 2.461964; // use known k
    double sum = 0;
    double n = input.size();
    for(size_t i = 0; i < input.size(); ++i) {
        sum += input[i];
    }
    double sigma = sum / (k * n);
    
    GammaParameters params;
    params.shape = k;
    params.rate = (1.0 / sigma);
    return params; 
}

double DurationModel::gamma_fit_calculate_s(const std::vector<double>& input)
{
    double sum_1 = 0;
    double sum_2 = 0;
    double n = input.size();
    for(size_t i = 0; i < input.size(); ++i) {
        sum_1 += input[i];
        assert(input[i] > 0.0f);
        sum_2 += log(input[i]);
    }

    return log(sum_1 / n) - sum_2 / n;
}
