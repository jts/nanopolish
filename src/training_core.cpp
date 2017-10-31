#include "training_core.hpp"
#include "nanopolish_emissions.h"
#include "logsumset.hpp"
#include "logger.hpp"

using std::string;
using std::vector;
using std::multiset;
using std::endl;

const bool use_multiset_logsum = 
#ifndef USE_MULTISET_LOGSUM
    false;
#else
    true;
#endif

ParamMixture train_gaussian_mixture(const vector< StateTrainingData >& data, const ParamMixture& input_mixture)
{
    size_t n_components = input_mixture.params.size();
    size_t n_data = data.size();
    float log_n_data = std::log(n_data);
    assert(input_mixture.log_weights.size() == n_components);
    ParamMixture curr_mixture = input_mixture;

    for(size_t iteration = 0; iteration < 10; ++iteration) {
        ParamMixture new_mixture = curr_mixture;

        // compute log_pdfs
        //
        //   pdf[i][j] := gauss(mu_j, sigma_j * read_var_i, level_mean_i)
        //
        vector< vector< float > > log_pdf(n_data);
        for(size_t i = 0; i < n_data; ++i) {
            log_pdf[i].resize(n_components);
            for(size_t j = 0; j < n_components; ++j) {
                // We need to scale the mixture component parameters by the per-read var factor
                PoreModelStateParams scaled_state = curr_mixture.params[j];
                scaled_state.level_stdv *= data[i].scaled_read_var;
                scaled_state.level_log_stdv += data[i].log_scaled_read_var;
                log_pdf[i][j] = log_normal_pdf(data[i].level_mean, scaled_state);
                assert(not std::isnan(log_pdf[i][j]));
                LOG("training_core", debug1)
                    << "pdf " << i << " " << j << " "
                    << std::scientific << std::exp(log_pdf[i][j]) << " ("
                    << std::fixed << std::setprecision(2) << log_pdf[i][j] << ")" << endl;
            }
        }

        // compute responsibilities
        //
        //   resp[i][j] := ( w_j * pdf[i][j] ) / sum_k ( w_k * pdf[i][k] )
        //
        vector< vector< float > > log_resp(n_data);
        for(size_t i = 0; i < n_data; ++i) {
            log_resp[i].resize(n_components);
            logsumset< float > denom_terms(use_multiset_logsum);
            for(size_t j = 0; j < n_components; ++j) {
                float v = curr_mixture.log_weights[j] + log_pdf[i][j];
                log_resp[i][j] = v;
                denom_terms.add(v);
                LOG("training_core", debug1)
                    << "resp_numer " << i << " " << j << " "
                    << std::scientific << std::exp(v) << " ("
                    << std::fixed << std::setprecision(2) << v << ")" << endl;
            }
            float log_denom = denom_terms.val();
            LOG("training_core", debug1) << "resp_denom " << i << " "
                << std::scientific << log_denom << " " << std::exp(log_denom) << std::endl;
            for(size_t j = 0; j < n_components; ++j) {
                log_resp[i][j] -= log_denom;
                LOG("training_core", debug1)
                    << "resp " << i << " " << j << " "
                    << std::fixed << std::setprecision(5) << std::exp(log_resp[i][j]) << " ("
                    << std::fixed << std::setprecision(2) << log_resp[i][j] << ")" << endl;
            }
        }

        // update weights
        //
        //   w'[j] := sum_i resp[i][j] / n_data
        //
        for (size_t j = 0; j < n_components; ++j) {
            logsumset< float > numer_terms(use_multiset_logsum);
            for (size_t i = 0; i < n_data; ++i) {
                numer_terms.add(log_resp[i][j]);
            }
            float log_numer = numer_terms.val();
            new_mixture.log_weights[j] = log_numer - log_n_data;
        }

        // update means
        //
        //   mu_j := sum_i ( resp[i][j] * level_mean_i ) / sum_i resp[i][j]
        //         = sum_i ( resp[i][j] * level_mean_i ) / ( w'[j] * n_data )
        //
        vector< float > new_log_mean(2);
        for (size_t j = 0; j < n_components; ++j) {
            logsumset< float > numer_terms(use_multiset_logsum);
            for (size_t i = 0; i < n_data; ++i) {
                numer_terms.add(log_resp[i][j] + data[i].log_level_mean);
            }
            float log_numer = numer_terms.val();
            new_log_mean[j] = log_numer - (log_n_data + new_mixture.log_weights[j]);
        }

        // update stdvs
        //
        //   var_j := sum_i ( resp[i][j] * ( ( level_mean_i - mu_j ) / scaled_read_var_i )^2 ) / sum_i resp[i][j]
        //          = sum_i ( resp[i][j] * ( ( level_mean_i - mu_j ) / scaled_read_var_i )^2 ) / ( w'[j] * n_data )
        //
        vector< float > new_log_var(2);
        for (size_t j = 0; j < n_components; ++j) {
            logsumset< float > numer_terms(use_multiset_logsum);
            for (size_t i = 0; i < n_data; ++i) {
                float v = std::abs(data[i].level_mean - std::exp(new_log_mean[j]));
                numer_terms.add(log_resp[i][j] + (not std::isnan(v) and v > 0? 2.0 * (std::log(v) - data[i].log_scaled_read_var) : 0.0));
            }
            float log_numer = numer_terms.val();
            new_log_var[j] = log_numer - (log_n_data + new_mixture.log_weights[j]);
        }

        for(size_t j = 0; j < n_components; ++j) {
            new_mixture.params[j].level_mean = std::exp(new_log_mean[j]);
            new_mixture.params[j].level_log_stdv = .5 * new_log_var[j];
            new_mixture.params[j].level_stdv = std::exp(new_mixture.params[j].level_log_stdv);
            LOG("training_core", debug)
                << "new_mixture " << iteration << " " << j << " "
                << std::fixed << std::setprecision(5) << std::exp(new_mixture.log_weights[j]) << " "
                << std::setprecision(3) << new_mixture.params[j].level_mean << " "
                << new_mixture.params[j].level_stdv << endl;
        }

        curr_mixture = new_mixture;
    }
    return curr_mixture;
}

ParamMixture train_invgaussian_mixture(const vector< StateTrainingData >& data, const ParamMixture& in_mixture)
{
    size_t n_components = in_mixture.params.size();
    assert(in_mixture.log_weights.size() == n_components);
    size_t n_data = data.size();
    auto crt_mixture = in_mixture;
    assert(false && "deprecated");
#if 0
    for (size_t j = 0; j < n_components; ++j) {
        LOG("training_core", debug)
            << "in_mixture " << j << " "
            << std::fixed << std::setprecision(5) << std::exp(in_mixture.log_weights[j]) << " "
            << std::setprecision(5) << in_mixture.params[j].sd_mean << endl;
    }

    // compute gaussian pdfs
    //
    //   pdf[i][j].first = gauss(mu_j, sigma_j * scaled_read_var_i, level_mean_i)
    //
    vector< vector< std::pair< float, float > > > log_pdf(n_data);
    for (size_t i = 0; i < n_data; ++i) {
        log_pdf[i].resize(n_components);
        for (size_t j = 0; j < n_components; ++j) {
            PoreModelStateParams scaled_state = in_mixture.params[j];
            scaled_state.level_stdv *= data[i].scaled_read_var;
            scaled_state.level_log_stdv += data[i].log_scaled_read_var;
            log_pdf[i][j].first = log_normal_pdf(data[i].level_mean, scaled_state);
            assert(not std::isnan(log_pdf[i][j].first));
            LOG("training_core", debug1)
                << "gauss_pdf " << i << " " << j << " "
                << std::scientific << std::exp(log_pdf[i][j].first) << " ("
                << std::fixed << std::setprecision(2) << log_pdf[i][j].first << ")" << endl;
        }
    }

    // compute gaussian weights
    //
    //   g_weights[i][j] := ( w_j * pdf[i][j].first ) / sum_k ( w_k * pdf[i][k].first )
    //
    vector< vector< float > > log_g_weights(n_data);
    for (size_t i = 0; i < n_data; ++i) {
        log_g_weights[i].resize(n_components);
        logsumset< float > denom_terms(use_multiset_logsum);
        for (size_t j = 0; j < n_components; ++j) {
            float v = in_mixture.log_weights[j] + log_pdf[i][j].first;
            log_g_weights[i][j] = v;
            denom_terms.add(v);
        }
        float log_denom = denom_terms.val();
        for (size_t j = 0; j < n_components; ++j) {
            log_g_weights[i][j] -= log_denom;
            LOG("training_core", debug1)
                << "g_weights " << i << " " << j << " "
                << std::fixed << std::setprecision(5) << std::exp(log_g_weights[i][j]) << " ("
                << std::fixed << std::setprecision(2) << log_g_weights[i][j] << ")" << endl;
        }
    }

    for (size_t iteration = 0; iteration < 10; ++iteration) {
        // compute inverse gaussian pdfs
        //
        //   pdf[i][j].second = invgauss(eta_j, lambda_j * ( read_var_sd_i / read_var_scale_i ), level_stdv_i)
        //
        for (size_t i = 0; i < n_data; ++i) {
            for (size_t j = 0; j < n_components; ++j) {
                PoreModelStateParams scaled_state = crt_mixture.params[j];
                scaled_state.sd_lambda *= data[i].read_var_sd / data[i].read_scale_sd;
                scaled_state.sd_log_lambda += data[i].log_read_var_sd - data[i].log_read_scale_sd;
                log_pdf[i][j].second = log_invgauss_pdf(data[i].level_stdv, data[i].log_level_stdv, scaled_state);
                assert(not std::isnan(log_pdf[i][j].second));
                LOG("training_core", debug1)
                    << "invgauss_pdf " << i << " " << j << " "
                    << std::scientific << std::exp(log_pdf[i][j].second) << " ("
                    << std::fixed << std::setprecision(2) << log_pdf[i][j].second << ")" << endl;
            }
        }
        // compute inverse gaussian weights (responsibilities)
        //
        //   ig_weights[i][j] := ( g_weights[i][j] * pdf[i][j].second ) / sum_k ( g_weights[i][k] * pdf[i][k].second )
        //
        vector< vector< float > > log_ig_weights(n_data);
        for (size_t i = 0; i < n_data; ++i) {
            log_ig_weights[i].resize(n_components);
            logsumset< float > denom_terms(use_multiset_logsum);
            for (size_t j = 0; j < n_components; ++j) {
                float v = log_g_weights[i][j] + log_pdf[i][j].second;
                log_ig_weights[i][j] = v;
                denom_terms.add(v);
            }
            float log_denom = denom_terms.val();
            for (size_t j = 0; j < n_components; ++j) {
                log_ig_weights[i][j] -= log_denom;
                LOG("training_core", debug1)
                    << "ig_weights " << i << " " << j << " "
                    << std::fixed << std::setprecision(5) << std::exp(log_ig_weights[i][j]) << " ("
                    << std::fixed << std::setprecision(2) << log_ig_weights[i][j] << ")" << endl;
            }
        }

        // update eta
        //
        //   eta_j := sum_i ( ig_weigts[i][j] * lambda'_ij * level_stdv_i ) / sum_i ( ig_weights[i][j] * lambda'_ij )
        //   lambda'_ij := lambda_j * ( read_var_sd_i / read_var_scale_i )
        //
        auto new_mixture = crt_mixture;
        for (size_t j = 0; j < n_components; ++j) {
            logsumset< float > numer_terms(use_multiset_logsum);
            logsumset< float > denom_terms(use_multiset_logsum);
            for (size_t i = 0; i < n_data; ++i) {
                float v = log_ig_weights[i][j] + crt_mixture.params[j].sd_log_lambda + (data[i].log_read_var_sd - data[i].log_read_scale_sd);
                numer_terms.add(v + data[i].log_level_stdv);
                denom_terms.add(v);
            }
            float log_numer = numer_terms.val();
            float log_denom = denom_terms.val();
            new_mixture.params[j].sd_mean = std::exp(log_numer - log_denom);
            new_mixture.params[j].update_sd_stdv();
            new_mixture.params[j].update_logs();
            LOG("training_core", debug)
                << "new_mixture " << iteration << " " << j << " "
                << std::fixed << std::setprecision(5) << std::exp(new_mixture.log_weights[j]) << " "
                << std::setprecision(5) << new_mixture.params[j].sd_mean << endl;
        }
        std::swap(crt_mixture, new_mixture);
    } // for iteration
#endif
    return crt_mixture;
} // train_ig_mixture
