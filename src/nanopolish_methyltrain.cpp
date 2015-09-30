//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_methyltrain -- train a methylation model
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include <cstddef>
#include "htslib/faidx.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

#include "../eigen/Eigen/Dense"

//
// Structs
//

// The state training data comes in two different
// sizes Full and Minimal. The model training functions
// only actually need the Minimal data but for exploration
// the Full data is useful so left as an option.
struct FullStateTrainingData
{
    //
    // Functions
    //
    FullStateTrainingData(const SquiggleRead& sr,
                          const EventAlignment& ea,
                          uint32_t rank,
                          const std::string& prev_kmer,
                          const std::string& next_kmer)
    {
        // scale the observation to the expected pore model
        this->level_mean = sr.get_fully_scaled_level(ea.event_idx, ea.strand_idx);
        //this->event_stdv = sr.events[strand_idx][ea.event_idx].stdv / sr.pore_model[strand_idx].scale_sd;
        this->level_stdv = 0;
        this->duration = sr.events[ea.strand_idx][ea.event_idx].duration;
        
        this->read_var = (float)sr.pore_model[ea.strand_idx].var;
        this->ref_position = ea.ref_position;
        this->ref_strand = ea.rc;
        
        GaussianParameters model = sr.pore_model[ea.strand_idx].get_scaled_parameters(rank);
        this->z = (sr.get_drift_corrected_level(ea.event_idx, ea.strand_idx) -  model.mean ) / model.stdv;
        this->prev_kmer = prev_kmer;
        this->next_kmer = next_kmer;
    }

    static void write_header(FILE* fp)
    {
        fprintf(fp, "model\tmodel_kmer\tlevel_mean\tlevel_stdv\tduration\tref_pos\tref_strand\tz\tread_var\tprev_kmer\tnext_kmer\n");
    }

    void write_tsv(FILE* fp, const std::string& model_name, const std::string& kmer) const
    {
        fprintf(fp, "%s\t%s\t%.2lf\t%.2lf\t%.3lf\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\n", 
                    model_name.c_str(), 
                    kmer.c_str(), 
                    level_mean, 
                    level_stdv,
                    duration,
                    ref_position,
                    ref_strand,
                    z,
                    read_var,
                    prev_kmer.c_str(),
                    next_kmer.c_str());
    }

    //
    // Data
    //

    float level_mean;
    float level_stdv;
    float duration;
    float read_var;

    int ref_position;
    int ref_strand;

    float z;
    std::string prev_kmer;
    std::string next_kmer;
};

struct MinimalStateTrainingData
{
    //
    // Functions
    //
    MinimalStateTrainingData(const SquiggleRead& sr,
                             const EventAlignment& ea,
                             uint32_t rank,
                             const std::string& prev_kmer,
                             const std::string& next_kmer)
    {
        // scale the observation to the expected pore model
        this->level_mean = sr.get_fully_scaled_level(ea.event_idx, ea.strand_idx);
        this->read_var = (float)sr.pore_model[ea.strand_idx].var;
    }

    static void write_header(FILE* fp)
    {
        fprintf(fp, "model\tmodel_kmer\tlevel_mean\tread_var\n");
    }

    void write_tsv(FILE* fp, const std::string& model_name, const std::string& kmer) const
    {
        fprintf(fp, "%s\t%s\t%.2lf\t%.2lf\n",
                    model_name.c_str(), 
                    kmer.c_str(), 
                    level_mean, 
                    read_var);
    }

    //
    // Data
    //
    float level_mean;
    float read_var;
};

typedef MinimalStateTrainingData StateTrainingData;

struct StateSummary
{
    StateSummary() { num_matches = 0; num_skips = 0; num_stays = 0; }
    std::vector<StateTrainingData> events;

    int num_matches;
    int num_skips;
    int num_stays;
};

struct GaussianMixture
{
    std::vector<float> weights;
    std::vector<GaussianParameters> params;
};

//
Alphabet* mtrain_alphabet = &gMCpGAlphabet;

//
// Typedefs
//
typedef std::map<std::string, std::vector<StateSummary>> ModelTrainingMap;

//
// Getopt
//
#define SUBPROGRAM "methyltrain"

static const char *METHYLTRAIN_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *METHYLTRAIN_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Train a methylation model\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -m, --models-fofn=FILE               read the models to be trained from the FOFN\n"
"      --train-unmethylated             train unmethylated 5-mers instead of methylated\n"
"  -c  --calibrate                      recalibrate aligned reads to model before training\n"
"      --no-update-models               do not write out trained models\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"  -s, --out-suffix=STR                 name output files like model.out_suffix\n"
"      --progress                       print out a progress message\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static unsigned int calibrate=0;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string models_fofn;
    static std::string region;
    static std::string out_suffix = ".methyltrain";
    static bool write_models = true;
    static bool train_unmethylated = false;
    static int progress = 0;
    static int num_threads = 1;
    static int batch_size = 128;
}

static const char* shortopts = "r:b:g:t:m:vnc";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PROGRESS, OPT_NO_UPDATE_MODELS, OPT_TRAIN_UNMETHYLATED };

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "calibrate",          no_argument,       NULL, 'c' },
    { "reads",              required_argument, NULL, 'r' },
    { "bam",                required_argument, NULL, 'b' },
    { "genome",             required_argument, NULL, 'g' },
    { "window",             required_argument, NULL, 'w' },
    { "threads",            required_argument, NULL, 't' },
    { "models-fofn",        required_argument, NULL, 'm' },
    { "out-suffix",         required_argument, NULL, 's' },
    { "no-update-models",   no_argument,       NULL, OPT_NO_UPDATE_MODELS },
    { "train-unmethylated", no_argument,       NULL, OPT_TRAIN_UNMETHYLATED },
    { "progress",           no_argument,       NULL, OPT_PROGRESS },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

GaussianMixture train_gaussian_mixture(const std::vector<StateTrainingData>& data,
                                       const GaussianMixture& input_mixture)
{

    size_t n_components = input_mixture.params.size();
    size_t n_data = data.size();
    assert(input_mixture.weights.size() == n_components);
    GaussianMixture curr_mixture = input_mixture;   

    for(size_t iteration = 0; iteration < 10; ++iteration) {
        std::vector<double> mean_sum(n_components, 0.0f);
        std::vector<double> var_sum(n_components, 0.0f);
        
        GaussianMixture new_mixture = curr_mixture;
        for(size_t j = 0; j < n_components; ++j) {
            new_mixture.weights[j] = 0.0f;
        }

        std::vector<std::vector<double> > resp;

        for(size_t i = 0; i < n_data; ++i) {
            // Calculate the posterior probability that
            // data i came from each component of the mixture
            
            // P(data i | component j) P(component j)
            std::vector<double> t(n_components, 0.0f);
            double t_sum = -INFINITY;
            for(size_t j = 0; j < n_components; ++j) {
                t[j] = log_normal_pdf(data[i].level_mean, curr_mixture.params[j]) + log(curr_mixture.weights[j]);
                if(t[j] != -INFINITY && ! std::isnan(t[j])) {
                    t_sum = add_logs(t_sum, t[j]);
                }
            }

            // store P(component j | data i)
            for(size_t j = 0; j < n_components; ++j) {
                t[j] = exp(t[j] - t_sum);
                new_mixture.weights[j] += t[j];
            }
            resp.push_back(t);
        }
        
        for(size_t j = 0; j < n_components; ++j) {
            new_mixture.weights[j] /= n_data;
        }

        // Calculate mean
        for(size_t i = 0; i < n_data; ++i) {
            for(size_t j = 0; j < n_components; ++j) {
                double w_ij = resp[i][j];
                mean_sum[j] += w_ij * data[i].level_mean;
            }
        }
        
        std::vector<double> new_mean(2);
        for(size_t j = 0; j < n_components; ++j) {
            new_mean[j] = mean_sum[j] / (n_data * new_mixture.weights[j]);
        }

        // Calculate variance
        for(size_t i = 0; i < n_data; ++i) {
            for(size_t j = 0; j < n_components; ++j) {
                double w_ij = resp[i][j];
                var_sum[j] += w_ij * pow( (data[i].level_mean - new_mean[j]) / data[i].read_var, 2.0);
            }
        }
        
        std::vector<double> new_var(2);
        for(size_t j = 0; j < n_components; ++j) {
            new_var[j] = var_sum[j] / (n_data * new_mixture.weights[j]);
        }

        for(size_t j = 0; j < n_components; ++j) {
            new_mixture.params[j] = GaussianParameters(new_mean[j], sqrt(new_var[j]));
            //fprintf(stderr, "MIXTURE\t%zu\t%.2lf\t%.2lf\t%.2lf\n", j, curr_mixture.weights[j], curr_mixture.params[j].mean, curr_mixture.params[j].stdv);
        }

        curr_mixture = new_mixture;
    }
    return curr_mixture;
}

double model_score(SquiggleRead &sr,
                   const size_t strand_idx,
                   const faidx_t *fai, 
                   const std::vector<EventAlignment> &alignment_output,
                   const size_t events_per_segment)  
{
    double curr_score = 0;
    size_t nevents = 0;

    for(int align_start_idx = events_per_segment; 
               align_start_idx < (int)alignment_output.size() - (int)events_per_segment; 
               align_start_idx += events_per_segment) {

        const EventAlignment& align_start = alignment_output[align_start_idx];
        const EventAlignment& align_end = alignment_output[align_start_idx + events_per_segment];
        std::string contig = alignment_output.front().ref_name.c_str();

        // Set up event data
        HMMInputData data;
        data.read = &sr;
        data.anchor_index = -1; // unused
        data.strand = strand_idx;
        data.rc = alignment_output.front().rc;
        data.event_start_idx = align_start.event_idx;
        data.event_stop_idx = align_end.event_idx;
        data.event_stride = data.event_start_idx <= data.event_stop_idx ? 1 : -1;
        
        // Set up reference data
        int ref_start_pos = align_start.ref_position;
        int ref_end_pos = align_end.ref_position;
        int fetched_len = 0;

        assert(ref_end_pos >= ref_start_pos);

        // Extract the reference sequence for this region
        std::string ref_seq = get_reference_region_ts(fai, contig.c_str(), ref_start_pos, 
                                                      ref_end_pos, &fetched_len);

        if (fetched_len < 100)
            continue;

        const Alphabet *alphabet = sr.pore_model[strand_idx].pmalphabet;
    
        ref_seq = alphabet->disambiguate(ref_seq);
        HMMInputSequence sequence(ref_seq, alphabet->reverse_complement(ref_seq), alphabet);

        // Run HMM using current model
        curr_score += profile_hmm_score(sequence, data, 0);
        nevents += events_per_segment;
    }

    if (nevents == 0)
        return +1;
    else
        return curr_score/nevents;
}

// recalculate shift, scale, drift, scale_sd from an alignment and the read
void recalibrate_model(SquiggleRead &sr,
                       const int strand_idx,
                       const std::vector<EventAlignment> &alignment_output, 
                       bool scale_var) 
{
    std::vector<double> raw_events, times, level_means, level_stdvs;
    uint32_t k = sr.pore_model[strand_idx].k;

    //std::cout << "Previous pore model parameters: " << sr.pore_model[strand_idx].shift << ", " 
    //                                                << sr.pore_model[strand_idx].scale << ", " 
    //                                                << sr.pore_model[strand_idx].drift << ", " 
    //                                                << sr.pore_model[strand_idx].var   << std::endl;

    // extract necessary vectors from the read and the pore model; note do not want scaled values
    for ( const auto &ea : alignment_output ) {
        if(ea.hmm_state == 'M') {
            std::string model_kmer = ea.rc ? mtrain_alphabet->reverse_complement(ea.ref_kmer) : ea.ref_kmer;
            uint32_t rank = mtrain_alphabet->kmer_rank(model_kmer.c_str(), k);

            raw_events.push_back ( sr.get_uncorrected_level(ea.event_idx, strand_idx) );
            times.push_back      ( sr.get_time(ea.event_idx, strand_idx) );
            level_means.push_back( sr.pore_model[strand_idx].states[rank].level_mean );
            level_stdvs.push_back( sr.pore_model[strand_idx].states[rank].level_stdv );
        }
    }

    const int minNumEventsToRescale = 500;
    if (raw_events.size() < minNumEventsToRescale) 
        return;

    // Assemble linear system corresponding to weighted least squares problem
    // Can just directly call a weighted least squares solver, but there's enough
    // structure in our problem it's a little faster just to build the normal eqn
    // matrices ourselves
    Eigen::Matrix3d A;
    Eigen::Vector3d b;

    for (int i=0; i<3; i++) {
        b(i) = 0.;
        for (int j=0; j<3; j++)
            A(i,j) = 0.;
    }

    for (size_t i=0; i<raw_events.size(); i++) {
        double inv_var = 1./(level_stdvs[i]*level_stdvs[i]);
        double mu = level_means[i];
        double t  = times[i];
        double e  = raw_events[i];

        A(0,0) += inv_var;  A(0,1) += mu*inv_var;    A(0,2) += t*inv_var;
                            A(1,1) += mu*mu*inv_var; A(1,2) += mu*t*inv_var;
                                                     A(2,2) += t*t*inv_var;

        b(0) += e*inv_var;
        b(1) += mu*e*inv_var;
        b(2) += t*e*inv_var;
    }
    A(1,0) = A(0,1);
    A(2,0) = A(0,2);
    A(2,1) = A(1,2);

    // perform the linear solve
    Eigen::Vector3d x = A.fullPivLu().solve(b);

    double shift = x(0);
    double scale = x(1);
    double drift = x(2);

    sr.pore_model[strand_idx].shift = shift;
    sr.pore_model[strand_idx].scale = scale;
    sr.pore_model[strand_idx].drift = drift;

    if (scale_var) {
        double var = 0.;
        for (size_t i=0; i<raw_events.size(); i++) {
            double yi = (raw_events[i] - shift - scale*level_means[i] - drift*times[i]);
            var+= yi*yi/(level_stdvs[i]*level_stdvs[i]);
        }
        var /= raw_events.size();

        sr.pore_model[strand_idx].var   = var;
    }

    if (sr.pore_model[strand_idx].is_scaled)
        sr.pore_model[strand_idx].bake_gaussian_parameters();

    //std::cout << "Updated pore model parameters:  " << sr.pore_model[strand_idx].shift << ", " 
    //                                                << sr.pore_model[strand_idx].scale << ", " 
    //                                                << sr.pore_model[strand_idx].drift << ", " 
    //                                                << sr.pore_model[strand_idx].var   << std::endl;
}

// Realign the read in event space
void train_read(const ModelMap& model_map,
                const Fast5Map& name_map, 
                const faidx_t* fai, 
                const bam_hdr_t* hdr, 
                const bam1_t* record, 
                size_t read_idx,
                int region_start,
                int region_end,
                size_t round,
                ModelTrainingMap& training)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);
    std::string fast5_path = name_map.get_path(read_name);

    // load read
    SquiggleRead sr(read_name, fast5_path);

    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {
        
        // replace model with the training model
        std::string curr_model = sr.pore_model[strand_idx].name;
        auto model_iter = model_map.find(curr_model);

        if(model_iter != model_map.end()) {
            sr.pore_model[strand_idx].update_states(model_iter->second);
        } else {
            printf("Error: model %s not found\n", curr_model.c_str());
            assert(false && "Model not found");
        }
        
        // set k
        uint32_t k = sr.pore_model[strand_idx].k;

        // Align to the new model
        EventAlignmentParameters params;
        params.sr = &sr;
        params.fai = fai;
        params.hdr = hdr;
        params.record = record;
        params.strand_idx = strand_idx;
 
        params.alphabet = mtrain_alphabet;
        params.read_idx = read_idx;
        params.region_start = region_start;
        params.region_end = region_end;
        std::vector<EventAlignment> alignment_output = align_read_to_ref(params);
        if (alignment_output.size() == 0)
            return;

        // Update pore model based on alignment
        double orig_score = model_score(sr, strand_idx, fai, alignment_output, 500);
#pragma omp critical(print)
        std::cout << round << " " << curr_model << " " << read_idx << " Original " << orig_score << std::endl;

        if ( opt::calibrate ) {
            recalibrate_model(sr, strand_idx, alignment_output, round > 5);

            double rescaled_score = model_score(sr, strand_idx, fai, alignment_output, 500);
#pragma omp critical(print)
            {
                std::cout << round << " " << curr_model << " " << read_idx << " Rescaled " << rescaled_score << std::endl;
                std::cout << round << " " << curr_model << " " << read_idx << " Delta " << rescaled_score-orig_score << std::endl;
            }
        }
        // Update model observations
        {
            //emit_event_alignment_tsv(stdout, sr, params, alignment_output);

            // Get the training data for this model
            auto& emission_map = training[curr_model];
            for(size_t i = 0; i < alignment_output.size(); ++i) {
                const EventAlignment& ea = alignment_output[i];
                std::string model_kmer = ea.model_kmer;

                // Grab the previous/next model kmer
                // If the read is from the same strand as the reference
                // the next kmer comes from the next alignment_output (and vice-versa)
                // other the indices are swapped
                int next_stride = ea.rc ? -1 : 1;

                std::string prev_kmer = "";
                std::string next_kmer = "";

                if(i > 0 && i < alignment_output.size() - 1) {
                    assert(alignment_output[i + next_stride].event_idx - ea.event_idx == 1);
                    assert(alignment_output[i - next_stride].event_idx - ea.event_idx == -1);

                    // check for exactly one base of movement along the reference
                    if( abs(alignment_output[i + next_stride].ref_position - ea.ref_position) == 1) {
                        next_kmer = alignment_output[i + next_stride].model_kmer;
                    }

                    if( abs(alignment_output[i - next_stride].ref_position - ea.ref_position) == 1) {
                        prev_kmer = alignment_output[i - next_stride].model_kmer;
                    }
                }

                uint32_t rank = mtrain_alphabet->kmer_rank(model_kmer.c_str(), k);
                auto& kmer_summary = emission_map[rank];
                
                // Should we use this event for training?
                bool use_for_training = i > 5 && 
                                        i + 5 < alignment_output.size() &&
                                        alignment_output[i].hmm_state == 'M' &&
                                        alignment_output[i - 1].hmm_state == 'M' &&
                                        alignment_output[i + 1].hmm_state == 'M' &&
                                        prev_kmer != "" &&
                                        next_kmer != "";

                if(use_for_training) {
                    StateTrainingData std(sr, ea, rank, prev_kmer, next_kmer);
                    #pragma omp critical(kmer)
                    kmer_summary.events.push_back(std);
                }

                if(ea.hmm_state == 'M')  {
                    #pragma omp atomic
                    kmer_summary.num_matches += 1;
                } else if(ea.hmm_state == 'E') {
                    #pragma omp atomic
                    kmer_summary.num_stays += 1;
                }

            }
        }
    } // for strands
}

ModelMap read_models_fofn(const std::string& fofn_name)
{
    ModelMap out;
    std::ifstream fofn_reader(fofn_name);
    std::string model_filename;

    while(getline(fofn_reader, model_filename)) {
        printf("reading %s\n", model_filename.c_str());
        PoreModel p(model_filename, *mtrain_alphabet);
        assert(!p.name.empty());

        out[p.name] = p;
    }
    return out;
}

void parse_methyltrain_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'm': arg >> opt::models_fofn; break;
            case 's': arg >> opt::out_suffix; break;
            case 'v': opt::verbose++; break;
            case 'c': opt::calibrate = 1; break;
            case OPT_TRAIN_UNMETHYLATED: opt::train_unmethylated = true; break;
            case OPT_NO_UPDATE_MODELS: opt::write_models = false; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << METHYLTRAIN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << METHYLTRAIN_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if(argc - optind > 0) {
        opt::region = argv[optind++];
    }

    if (argc - optind > 0) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::num_threads <= 0) {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::num_threads << "\n";
        die = true;
    }

    if(opt::reads_file.empty()) {
        std::cerr << SUBPROGRAM ": a --reads file must be provided\n";
        die = true;
    }
    
    if(opt::genome_file.empty()) {
        std::cerr << SUBPROGRAM ": a --genome file must be provided\n";
        die = true;
    }

    if(opt::bam_file.empty()) {
        std::cerr << SUBPROGRAM ": a --bam file must be provided\n";
        die = true;
    }
    
    if(opt::models_fofn.empty()) {
        std::cerr << SUBPROGRAM ": a --models-fofn file must be provided\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << METHYLTRAIN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

ModelMap train_one_round(const ModelMap& models, const Fast5Map& name_map, size_t round)
{
    // Initialize the training summary stats for each kmer for each model
    ModelTrainingMap model_training_data;
    for(auto model_iter = models.begin(); model_iter != models.end(); model_iter++) {
        std::vector<StateSummary> summaries(model_iter->second.get_num_states()); // one per kmer in the model
        model_training_data[model_iter->first] = summaries;
    }

    // Open the BAM and iterate over reads

    // load bam file
    htsFile* bam_fh = sam_open(opt::bam_file.c_str(), "r");
    assert(bam_fh != NULL);

    // load bam index file
    std::string index_filename = opt::bam_file + ".bai";
    hts_idx_t* bam_idx = bam_index_load(index_filename.c_str());
    assert(bam_idx != NULL);

    // read the bam header
    bam_hdr_t* hdr = sam_hdr_read(bam_fh);
    
    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    hts_itr_t* itr;

    // If processing a region of the genome, only emit events aligned to this window
    int clip_start = -1;
    int clip_end = -1;

    if(opt::region.empty()) {
        // TODO: is this valid?
        itr = sam_itr_queryi(bam_idx, HTS_IDX_START, 0, 0);
    } else {

        fprintf(stderr, "Region: %s\n", opt::region.c_str());
        itr = sam_itr_querys(bam_idx, hdr, opt::region.c_str());
        hts_parse_reg(opt::region.c_str(), &clip_start, &clip_end);
    }

#ifndef H5_HAVE_THREADSAFE
    if(opt::num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif

    // Initialize iteration
    std::vector<bam1_t*> records(opt::batch_size, NULL);
    for(size_t i = 0; i < records.size(); ++i) {
        records[i] = bam_init1();
    }

    int result;
    size_t num_reads_realigned = 0;
    size_t num_records_buffered = 0;
    Progress progress("[methyltrain]");

    do {
        assert(num_records_buffered < records.size());
        
        // read a record into the next slot in the buffer
        result = sam_itr_next(bam_fh, itr, records[num_records_buffered]);
        num_records_buffered += result >= 0;

        // realign if we've hit the max buffer size or reached the end of file
        if(num_records_buffered == records.size() || result < 0) {
            #pragma omp parallel for            
            for(size_t i = 0; i < num_records_buffered; ++i) {
                bam1_t* record = records[i];
                size_t read_idx = num_reads_realigned + i;
                if( (record->core.flag & BAM_FUNMAP) == 0) {
                    train_read(models, name_map, fai, hdr, record, read_idx, clip_start, clip_end, round, model_training_data);
                }
            }

            num_reads_realigned += num_records_buffered;
            num_records_buffered = 0;
        }

        if(opt::progress) {
            fprintf(stderr, "Realigned %zu reads in %.1lfs\r", num_reads_realigned, progress.get_elapsed_seconds());
        }
    } while(result >= 0);
    
    assert(num_records_buffered == 0);
    progress.end();

    std::stringstream training_fn;
    training_fn << opt::bam_file << ".methyltrain.tsv";
    
    std::stringstream summary_fn;
    summary_fn << opt::bam_file << ".methyltrain.summary";

    FILE* summary_fp = fopen(summary_fn.str().c_str(), "w");
    FILE* training_fp = fopen(training_fn.str().c_str(), "w");

    // training header
    StateTrainingData::write_header(training_fp);
    
    // Process the training results
    ModelMap trained_models;
    
    for(auto model_training_iter = model_training_data.begin(); 
             model_training_iter != model_training_data.end(); model_training_iter++) {
        
        // Initialize trained model from input model
        auto model_iter = models.find(model_training_iter->first);
        assert(model_iter != models.end());

        std::string model_name = model_training_iter->first;
        std::string model_short_name = "";

        if(model_name == "r7.3_template_median68pA.model") {
            model_short_name = "t";
        } else if(model_name == "r7.3_complement_median68pA_pop1.model") {
            model_short_name = "c.p1";
        } else if(model_name == "r7.3_complement_median68pA_pop2.model") {
            model_short_name = "c.p2";
        } else if(model_name == "r7.3_e6_70bps_6mer_template_median68pA.model") {
            model_short_name = "t.006";
        } else if(model_name == "r7.3_e6_70bps_6mer_complement_median68pA_pop1.model") {
            model_short_name = "c.p1.006";
        } else if(model_name == "r7.3_e6_70bps_6mer_complement_median68pA_pop2.model") {
            model_short_name = "c.p2.006";
        } else {
            printf("Unknown model: %s\n", model_name.c_str());
            assert(false);
        }

        trained_models[model_training_iter->first] = model_iter->second;
        PoreModel& new_pm = trained_models[model_training_iter->first];

        // Update means for each kmer
        uint32_t k = new_pm.k;
        std::string kmer(k, 'A');
        const std::vector<StateSummary>& summaries = model_training_iter->second;
        for(size_t ki = 0; ki < summaries.size(); ++ki) {

            // write a training file
            for(size_t ei = 0; ei < summaries[ki].events.size(); ++ei) {
                summaries[ki].events[ei].write_tsv(training_fp, model_short_name, kmer);
            }

            // write to the summary file
            fprintf(summary_fp, "%s\t%s\t%d\t%d\t%d\n", model_short_name.c_str(), kmer.c_str(), summaries[ki].num_matches, summaries[ki].num_skips, summaries[ki].num_stays);

            GaussianMixture mixture;

            // train a mixture model where a minority of k-mers aren't methylated
            
            // unmethylated component
            double um_rate = 0.05f;
            std::string um_kmer = gMCpGAlphabet.unmethylate(kmer);
            size_t um_ki = gMCpGAlphabet.kmer_rank(um_kmer.c_str(), k);
            GaussianParameters um_params(model_iter->second.get_parameters(um_ki).level_mean, 
                                           model_iter->second.get_parameters(um_ki).level_stdv);

            mixture.weights.push_back(um_rate);
            mixture.params.push_back(um_params);

            GaussianParameters m_params(model_iter->second.get_parameters(ki).level_mean, 
                                           model_iter->second.get_parameters(ki).level_stdv);

            mixture.weights.push_back(1 - um_rate);
            mixture.params.push_back(m_params);
 
            if(opt::verbose > 1) {
                           
                fprintf(stderr, "INIT__MIX %s\t%s\t[%.2lf %.2lf %.2lf]\t[%.2lf %.2lf %.2lf]\n", model_training_iter->first.c_str(), kmer.c_str(), 
                                                                  mixture.weights[0], mixture.params[0].mean, mixture.params[0].stdv,
                                                                  mixture.weights[1], mixture.params[1].mean, mixture.params[1].stdv);
            }

            GaussianMixture trained_mixture = train_gaussian_mixture(summaries[ki].events, mixture);
            
            if(opt::verbose > 1) {
                fprintf(stderr, "TRAIN_MIX %s\t%s\t[%.2lf %.2lf %.2lf]\t[%.2lf %.2lf %.2lf]\n", model_training_iter->first.c_str(), kmer.c_str(), 
                                                                  trained_mixture.weights[0], trained_mixture.params[0].mean, trained_mixture.params[0].stdv,
                                                                  trained_mixture.weights[1], trained_mixture.params[1].mean, trained_mixture.params[1].stdv);
            }
                
            bool is_m_kmer = kmer.find('M') != std::string::npos;
            bool update_kmer = (is_m_kmer == !opt::train_unmethylated); 
            if(update_kmer && summaries[ki].events.size() > 100) {
                new_pm.states[ki].level_mean = trained_mixture.params[1].mean;
                new_pm.states[ki].level_stdv = trained_mixture.params[1].stdv;
            }

            /*
            if(kmer.find("CG") != std::string::npos) {
                float mu_prime = summaries[ki].mean_sum / summaries[ki].n;
                float var_prime = summaries[ki].var_sum / summaries[ki].n;
                new_pm[ki].level_mean = mu_prime;
                new_pm[ki].level_stdv = sqrt(var_prime);
                fprintf(stderr, "%s %s %.2lf %.2lf\n", model_training_iter->first.c_str(), kmer.c_str(), new_pm[ki].level_mean, new_pm[ki].level_stdv);
            }
            */
            mtrain_alphabet->lexicographic_next(kmer);
        }
    }


    // cleanup records
    for(size_t i = 0; i < records.size(); ++i) {
        bam_destroy1(records[i]);
    }

    // cleanup
    sam_itr_destroy(itr);
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
    fclose(training_fp);
    fclose(summary_fp);
    return trained_models;
}

void write_models(ModelMap& models)
{
    // Write the model
    for(auto model_iter = models.begin(); 
             model_iter != models.end(); model_iter++) {

        assert(!model_iter->second.model_filename.empty());
        std::string outname   =  model_iter->second.model_filename + opt::out_suffix;
        std::string modelname =  model_iter->first + (!opt::train_unmethylated ? opt::out_suffix : "");
        models[model_iter->first].write( outname, modelname );
        models[model_iter->first].write( outname, *mtrain_alphabet, modelname );
    }
}

int methyltrain_main(int argc, char** argv)
{
    parse_methyltrain_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    Fast5Map name_map(opt::reads_file);
    ModelMap models = read_models_fofn(opt::models_fofn);
    
    const size_t TRAINING_ROUNDS = 10;
    for(size_t round = 0; round < TRAINING_ROUNDS; round++) {
        fprintf(stderr, "Starting round %zu\n", round);
        ModelMap trained_models = train_one_round(models, name_map, round);
        if(opt::write_models) {
            write_models(trained_models);
        }
        models = trained_models;
    }
    return EXIT_SUCCESS;
}

