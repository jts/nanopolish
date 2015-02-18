// TODO: Boilerplate

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "hmmcons_khmm_parameters.h"

void khmm_parameters_initialize(KHMMParameters& parameters)
{
    parameters.skip_bin_width = 0.5;
    parameters.skip_probabilities.resize(40);

    // These default values are learned from a set of e.coli reads
    // trained on a de novo assembly
    parameters.self_transition = 0.2f;

    parameters.skip_probabilities[0] = 0.51268137;
    parameters.skip_probabilities[1] = 0.47243219;
    parameters.skip_probabilities[2] = 0.42888741;
    parameters.skip_probabilities[3] = 0.34932588;
    parameters.skip_probabilities[4] = 0.27427068;
    parameters.skip_probabilities[5] = 0.22297225;
    parameters.skip_probabilities[6] = 0.17585147;
    parameters.skip_probabilities[7] = 0.14705882;
    parameters.skip_probabilities[8] = 0.12183525;
    parameters.skip_probabilities[9] = 0.11344997;
    parameters.skip_probabilities[10] = 0.10069393;
    parameters.skip_probabilities[11] = 0.09153005;
    parameters.skip_probabilities[12] = 0.08765206;
    parameters.skip_probabilities[13] = 0.08491435;
    parameters.skip_probabilities[14] = 0.08272553;
    parameters.skip_probabilities[15] = 0.07747396;
    parameters.skip_probabilities[16] = 0.08439116;
    parameters.skip_probabilities[17] = 0.07819045;
    parameters.skip_probabilities[18] = 0.07337461;
    parameters.skip_probabilities[19] = 0.07020490;
    parameters.skip_probabilities[20] = 0.06869961;
    parameters.skip_probabilities[21] = 0.06576609;
    parameters.skip_probabilities[22] = 0.06923376;
    parameters.skip_probabilities[23] = 0.06239092;
    parameters.skip_probabilities[24] = 0.06586513;
    parameters.skip_probabilities[25] = 0.07372986;
    parameters.skip_probabilities[26] = 0.07050360;
    parameters.skip_probabilities[27] = 0.07228916;
    parameters.skip_probabilities[28] = 0.05855856;
    parameters.skip_probabilities[29] = 0.06842737;
    parameters.skip_probabilities[30] = 0.06145251;
    parameters.skip_probabilities[31] = 0.07352941;
    parameters.skip_probabilities[32] = 0.06278027;
    parameters.skip_probabilities[33] = 0.05932203;
    parameters.skip_probabilities[34] = 0.09708738;
    parameters.skip_probabilities[35] = 0.08290155;
    parameters.skip_probabilities[36] = 0.07692308;
    parameters.skip_probabilities[37] = 0.06896552;
    parameters.skip_probabilities[38] = 0.03448276;
    parameters.skip_probabilities[39] = 0.02985075;

    // initialize training data
    parameters.training_data.n_matches = 0;
    parameters.training_data.n_merges = 0;
    parameters.training_data.n_skips = 0;

    parameters.fit_quality = 0.0f;
}

inline size_t get_bin(const KHMMParameters& parameters, double k_level1, double k_level2)
{
    assert(!parameters.skip_probabilities.empty());

    double d = fabs(k_level1 - k_level2);
    size_t bin = d / parameters.skip_bin_width;

    // clamp out-of-range to last value
    bin = bin >= parameters.skip_probabilities.size() ? parameters.skip_probabilities.size() - 1 : bin;
    return bin;
}

double get_skip_probability(const KHMMParameters& parameters, double k_level1, double k_level2)
{
    size_t bin = get_bin(parameters, k_level1, k_level2);
    assert(bin < parameters.skip_probabilities.size());
    return parameters.skip_probabilities[bin];
}

void khmm_parameters_train(KHMMParameters& parameters)
{
    TrainingData& td = parameters.training_data;
    printf("train -- M: %d E: %d K: %d\n", td.n_matches, td.n_merges, td.n_skips);

    size_t sum = td.n_matches + td.n_merges + td.n_skips;
    if(sum > 0)
        parameters.self_transition = (double)td.n_merges / sum;
    printf("SKIPLEARN -- self %.3lf\n", parameters.self_transition);

    // Initialize observations with pseudocounts from the current model
    size_t num_bins = parameters.skip_probabilities.size();
    uint32_t pseudocount = 100;
    std::vector<double> total_observations(num_bins, 0.0f);
    std::vector<double> skip_observations(num_bins, 0.0f);

    for(size_t bin = 0; bin < num_bins; bin++) {
        skip_observations[bin] = parameters.skip_probabilities[bin] * pseudocount;
        total_observations[bin] = pseudocount;
    }

    // Iterate over the observed data and update counts
    for(size_t oi = 0; oi < td.transitions.size(); ++oi) {
        const TransitionObservation& to = td.transitions[oi];
        bool is_skip = to.state == 'K';
        size_t bin = get_bin(parameters, to.level_1, to.level_2);
        skip_observations[bin] += is_skip;
        total_observations[bin] += 1;
    }

    // Update probabilities
    for(size_t bin = 0; bin < num_bins; bin++) {
        parameters.skip_probabilities[bin] = skip_observations[bin] / total_observations[bin];
        printf("SKIPLEARN -- bin[%zu] %.3lf %.3lf %.3lf\n", bin, skip_observations[bin], total_observations[bin], parameters.skip_probabilities[bin]);
    }

    // Calculate how well the emissions fit the expected model
    parameters.fit_quality = 0.0f;
    for(size_t ei = 0; ei < td.emissions_for_matches.size(); ++ei) {
        parameters.fit_quality += td.emissions_for_matches[ei];
    }

    if(td.emissions_for_matches.size() > 0)
        parameters.fit_quality /= td.emissions_for_matches.size();
    else
        parameters.fit_quality = -INFINITY;
    printf("MODELFIT: %.2lf\n", parameters.fit_quality);
}
