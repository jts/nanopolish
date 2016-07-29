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
    probability_by_index.resize(100, 0.0f);
    probability_by_index[0] = 0.04966;
    probability_by_index[4] = 0.07400;
    probability_by_index[6] = 0.05673;
    probability_by_index[7] = 0.05538;
    probability_by_index[8] = 0.10918;
    probability_by_index[9] = 0.05356;
    probability_by_index[10] = 0.05009;
    probability_by_index[11] = 0.00152;
    probability_by_index[12] = 0.04834;
    probability_by_index[13] = 0.05119;
    probability_by_index[14] = 0.06148;
    probability_by_index[15] = 0.02331;
    probability_by_index[16] = 0.03679;
    probability_by_index[17] = 0.02491;
    probability_by_index[18] = 0.02903;
    probability_by_index[19] = 0.02655;
    probability_by_index[20] = 0.02332;
    probability_by_index[21] = 0.02150;
    probability_by_index[22] = 0.00690;
    probability_by_index[23] = 0.01939;
    probability_by_index[24] = 0.01670;
    probability_by_index[25] = 0.01432;
    probability_by_index[26] = 0.01339;
    probability_by_index[27] = 0.01196;
    probability_by_index[28] = 0.00624;
    probability_by_index[29] = 0.01008;
    probability_by_index[30] = 0.00898;
    probability_by_index[31] = 0.00759;
    probability_by_index[32] = 0.00757;
    probability_by_index[33] = 0.00652;
    probability_by_index[34] = 0.00323;
    probability_by_index[35] = 0.00581;
    probability_by_index[36] = 0.00486;
    probability_by_index[37] = 0.00472;
    probability_by_index[38] = 0.00424;
    probability_by_index[39] = 0.00232;
    probability_by_index[40] = 0.00366;
    probability_by_index[41] = 0.00336;
    probability_by_index[42] = 0.00290;
    probability_by_index[43] = 0.00281;
    probability_by_index[44] = 0.00257;
    probability_by_index[45] = 0.00141;
    probability_by_index[46] = 0.00219;
    probability_by_index[47] = 0.00205;
    probability_by_index[48] = 0.00173;
    probability_by_index[49] = 0.00170;
    probability_by_index[50] = 0.00158;
    probability_by_index[51] = 0.00095;
    probability_by_index[52] = 0.00140;
    probability_by_index[53] = 0.00119;
    probability_by_index[54] = 0.00118;
    probability_by_index[55] = 0.00106;
    probability_by_index[56] = 0.00065;
    probability_by_index[57] = 0.00094;
    probability_by_index[58] = 0.00091;
    probability_by_index[59] = 0.00077;
    probability_by_index[60] = 0.00077;
    probability_by_index[61] = 0.00073;
    probability_by_index[62] = 0.00045;
    probability_by_index[63] = 0.00068;
    probability_by_index[64] = 0.00061;
    probability_by_index[65] = 0.00053;
    probability_by_index[66] = 0.00050;
    probability_by_index[67] = 0.00050;
    probability_by_index[68] = 0.00031;
    probability_by_index[69] = 0.00042;
    probability_by_index[70] = 0.00041;
    probability_by_index[71] = 0.00034;
    probability_by_index[72] = 0.00035;
    probability_by_index[73] = 0.00023;
    probability_by_index[74] = 0.00033;
    probability_by_index[75] = 0.00031;
    probability_by_index[76] = 0.00026;
    probability_by_index[77] = 0.00027;
    probability_by_index[78] = 0.00027;
    probability_by_index[79] = 0.00018;
    probability_by_index[80] = 0.00022;
    probability_by_index[81] = 0.00021;
    probability_by_index[82] = 0.00018;
    probability_by_index[83] = 0.00019;
    probability_by_index[84] = 0.00017;
    probability_by_index[85] = 0.00012;
    probability_by_index[86] = 0.00017;
    probability_by_index[87] = 0.00014;
    probability_by_index[88] = 0.00014;
    probability_by_index[89] = 0.00014;
    probability_by_index[90] = 0.00010;
    probability_by_index[91] = 0.00013;
    probability_by_index[92] = 0.00013;
    probability_by_index[93] = 0.00011;
    probability_by_index[94] = 0.00011;
    probability_by_index[95] = 0.00010;
    probability_by_index[96] = 0.00008;
    probability_by_index[97] = 0.00009;
    probability_by_index[98] = 0.00009;
    probability_by_index[99] = 0.00251;
}

void DurationModel::add_training_from_sequence(const std::string& sequence, 
                                               const HMMInputData& data,
                                               const uint32_t alignment_flags)
{
    DurationModel& instance = getInstance();
    std::vector<double> duration_by_kmer_position = generate_aligned_durations(sequence, data, alignment_flags);
    
    #pragma omp critical
    for(const auto& duration : duration_by_kmer_position) {
        instance.training_data.push_back(duration);
    }
}

std::vector<double> DurationModel::generate_aligned_durations(const std::string& sequence,
                                                              const HMMInputData& data,
                                                              const uint32_t alignment_flags)
{
    size_t k = data.read->pore_model[0].k;
    size_t num_kmers = sequence.size() - k + 1;
    // initialize the vector of durations
    std::vector<double> duration_by_kmer_position(num_kmers, 0.0);
    std::vector<HMMAlignmentState> alignment = profile_hmm_align(sequence, data, alignment_flags);
    for(size_t ai = 0; ai < alignment.size(); ai++) {

        if(alignment[ai].state == 'M') {
            duration_by_kmer_position[alignment[ai].kmer_idx] += data.read->get_duration(alignment[ai].event_idx, data.strand);
        }
    }
    return duration_by_kmer_position;
}

/*
void DurationModel::generate_aligned_durations(const std::string& sequence,
                                               const HMMInputData& data)
{
    std::vector<HMMAlignmentState> alignment = profile_hmm_align(sequence, data);

    // The HMM can skip k-mers by not aligning an event to them
    // We want to avoid having zero-duration k-mers in the model
    // so we distribute the durations of the events assigned to the k-mers
    // flanking the gap over the whole region.
    size_t prev_kmer_idx = -1;
    for(size_t ai = 0; ai < alignment.size(); ai++) {

        size_t curr_kmer_idx = alignment[ai].kmer_idx;

        // Get the index of the next kmer aligned to
        size_t end_ai = ai + 1;
        size_t next_kmer_idx = -1;
        while(end_ai < alignment.size() && alignment[end_ai].kmer_idx == curr_kmer_idx) {

        }
    }
}
*/

void DurationModel::print_training_data(FILE* fp)
{
    DurationModel& instance = getInstance();
    for(size_t i = 0; i < instance.training_data.size(); ++i) {
        fprintf(fp, "%.5lf\n", instance.training_data[i]);
    }
}

void DurationModel::print_compiled_model(FILE* fp)
{
    
    DurationModel& instance = getInstance();
    
    std::map<int, int> count_by_indexed_duration;
    int total = 0;
    for(size_t i = 0; i < instance.training_data.size(); ++i) {
        double duration = instance.training_data[i];
        int index = get_index(duration);
        count_by_indexed_duration[index] += 1;
        total += 1;
    }

    for(auto& pair : count_by_indexed_duration) {
        fprintf(fp, "%d %.5lf %d %.5lf\n", pair.first, pair.first * MIN_DURATION, pair.second, pair.second / (double)total);
    }
    
    for(auto& pair : count_by_indexed_duration) {
        fprintf(fp, "probability_by_index[%d] = %.5lf;\n", pair.first, pair.second / (double)total);
    }
    
}

double DurationModel::score_sequence(const std::string& sequence, const HMMInputData& data, const uint32_t alignment_flags)
{
    DurationModel& instance = getInstance();
    std::vector<double> duration_by_kmer_position = generate_aligned_durations(sequence, data, alignment_flags);
    
    double lp = 0.0;
    for(const auto& duration : duration_by_kmer_position) {
        lp += log(instance.probability_by_index[get_index(duration)]);
    }
    return lp;
}

std::vector<double> DurationModel::score_vector_sequence(const std::string& sequence, const HMMInputData& data, const uint32_t alignment_flags)
{
    DurationModel& instance = getInstance();
    std::vector<double> duration_by_kmer_position = generate_aligned_durations(sequence, data, alignment_flags);
    std::vector<double> out;
       
    for(const auto& duration : duration_by_kmer_position) {
        out.push_back(log(instance.probability_by_index[get_index(duration)]));
    }
    return out;
}

double DurationModel::log_gamma_sum(double x, double n)
{
    const static double a = 2.461964; // shape
    const static double b = 587.2858; // rate
    assert(x >= 0.0);

    double na = n * a;
    return (na * log(b)) - lgamma(na) + (na - 1) * log(x) - b * x;
}
