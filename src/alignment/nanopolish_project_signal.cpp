//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_project_signal -- project signal data onto a reference genome
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
#include <sstream>
#include <iomanip>
#include <set>
#include <omp.h>
#include <getopt.h>
#include <math.h>
#include <iterator>
#include <fstream>
#include "htslib/faidx.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_read_db.h"
#include "nanopolish_hmm_input_sequence.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_bam_processor.h"
#include "nanopolish_bam_utils.h"
#include "nanopolish_raw_loader.h"
#include "fs_support.hpp"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

//
using namespace std::placeholders;

//
// Getopt
//
#define SUBPROGRAM "project-signal"

static const char *PROJECT_SIGNAL_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2023 Ontario Institute for Cancer Research\n";

static const char *PROJECT_SIGNAL_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --directory fast5s/ --kmer-model levels.txt --bam alignments.bam --genome genome.fa\n"
"Align nanopore events to reference k-mers\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"  -m, --kmer-model-file                read expected current levels from FILE\n"
"  -d, --directory                      path to the directory containing the raw ONT signal files. This path will be searched recursively so can contain subdirs\n"
"  -q, --min-mapping-quality=NUM        only use reads with mapping quality at least NUM (default: 0)\n"
//"      --progress                       print out a progress message\n"
"      --samples                        write the raw samples for the event to the tsv output\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string bam_file;
    static std::string genome_file;
    static std::string region;
    static std::string summary_file;
    static std::string kmer_model_file;
    static std::string fast5_directory;
    static int progress = 0;
    static int num_threads = 1;
    static int batch_size = 512;
    static int min_mapping_quality = 0;
    static int alignment_end_trim = 20;
    static bool write_samples = false;
}

static const char* shortopts = "b:g:t:q:m:d:vn";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PROGRESS, OPT_SAM, OPT_SUMMARY, OPT_SCALE_EVENTS, OPT_MODELS_FOFN, OPT_SAMPLES, OPT_SIGNAL_INDEX };

static const struct option longopts[] = {
    { "verbose",             no_argument,       NULL, 'v' },
    { "bam",                 required_argument, NULL, 'b' },
    { "genome",              required_argument, NULL, 'g' },
    { "threads",             required_argument, NULL, 't' },
    { "min-mapping-quality", required_argument, NULL, 'q' },
    { "kmer-model-file",     required_argument, NULL, 'm' },
    { "directory",           required_argument, NULL, 'd' },
    { "samples",             no_argument,       NULL, OPT_SAMPLES },
    { "progress",            no_argument,       NULL, OPT_PROGRESS },
    { "help",                no_argument,       NULL, OPT_HELP },
    { "version",             no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

std::vector<AlignedPair> sample_dp(std::vector<float> expected_levels,
                                   std::vector<float> scaled_samples,
                                   std::vector<int> start_samples,
                                   float skip_penalty,
                                   float trim_penalty)
{
    std::vector<AlignedPair> dummy_return;
    AdaBandedParameters params;
    EventBandedViterbi hmm_result;
    hmm_result.initialize(expected_levels, scaled_samples, start_samples, params);

    size_t n_events = hmm_result.get_num_events();
    size_t n_kmers = hmm_result.get_num_kmers();

    //float trim_penalty = -INFINITY;
    float step_penalty = 0.0f;
    float stay_penalty = 0.0f;
    //float skip_penalty = 0.0f;

    assert(hmm_result.is_initialized());

    // initialize first two bands as a special case

    // set origin cell
    hmm_result.set3_by_event_kmer(-1, -1, 0.0f, -INFINITY, -INFINITY);

    // fill in remaining bands
    for(size_t band_idx = 1; band_idx < hmm_result.get_num_bands() - 1; ++band_idx) {

        hmm_result.determine_band_origin(band_idx);

        // update start trim state for this band
        int start_trim_kmer_state = -1;
        int start_trim_offset = hmm_result.get_offset_for_kmer_in_band(band_idx, start_trim_kmer_state);
        if(hmm_result.is_offset_valid(start_trim_offset)) {
            int event_idx = hmm_result.get_event_at_band_offset(band_idx, start_trim_offset);
            float score_u = hmm_result.get_by_event_kmer(event_idx - 1, start_trim_kmer_state) + trim_penalty;
            hmm_result.set3(band_idx, start_trim_offset, -INFINITY, score_u, -INFINITY);
        }
 
        int min_offset, max_offset;
        hmm_result.get_offset_range_for_band(band_idx, min_offset, max_offset);

        for(int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = hmm_result.get_event_at_band_offset(band_idx, offset);
            int kmer_idx = hmm_result.get_kmer_at_band_offset(band_idx, offset);

            float diag = hmm_result.get_by_event_kmer(event_idx - 1, kmer_idx - 1);
            float up   = hmm_result.get_by_event_kmer(event_idx - 1, kmer_idx);
            float left = hmm_result.get_by_event_kmer(event_idx, kmer_idx - 1);
            
#ifdef VERIFY_MEMORY
            assert(event_idx >= 0 && event_idx < n_events);
            assert(kmer_idx >= 0 && kmer_idx < n_kmers);
#endif
            // the viterbi framework expects to maximize the log score,
            // so we multiple by -1 to turn the minimization problem into a maximization
            float score = -1.0 * fabsf(scaled_samples[event_idx] - expected_levels[kmer_idx]);
            float score_d = diag + step_penalty + score;
            float score_u = up + stay_penalty + score;
            float score_l = left + skip_penalty;
            hmm_result.set3_by_event_kmer(event_idx, kmer_idx, score_d, score_u, score_l);

#ifdef DEBUG_GENERIC
            fprintf(stderr, "[ada-gen-fill] bi: %d o: %d e: %d k: %d s: %.2lf rank: %zu emit: %.2lf\n", 
                band_idx, offset, event_idx, kmer_idx, hmm_result.get(band_idx, offset), kmer_rank, lp_emission);
            fprintf(stderr, "[ada-gen-fill]\tup: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left);
#endif
        }

        // if there is an end trim state in this band, set it here
        int end_trim_kmer_state = n_kmers;
        int offset = hmm_result.get_offset_for_kmer_in_band(band_idx, end_trim_kmer_state);
        if(hmm_result.is_offset_valid(offset)) {
            int event_idx = hmm_result.get_event_at_band_offset(band_idx, offset);
            float score_d = hmm_result.get_by_event_kmer(event_idx - 1, n_kmers - 1) + step_penalty;
            float score_u = hmm_result.get_by_event_kmer(event_idx - 1, end_trim_kmer_state) + trim_penalty;
            float score_l = hmm_result.get_by_event_kmer(event_idx, n_kmers - 1) + skip_penalty;
            hmm_result.set3(band_idx, offset, score_d, score_u, score_l);
#ifdef DEBUG_GENERIC
            fprintf(stderr, "[ada-gen-fill] set end trim %zu %d %d %d\n", band_idx, end_trim_kmer_state, event_idx, offset);
#endif
        }
    }

    // terminate
    int terminal_event_idx = n_events;
    int terminal_kmer_idx = n_kmers;

    float score_d = hmm_result.get_by_event_kmer(terminal_event_idx - 1, terminal_kmer_idx - 1);
    float score_u = hmm_result.get_by_event_kmer(terminal_event_idx - 1, terminal_kmer_idx);
    float score_l = hmm_result.get_by_event_kmer(terminal_event_idx, terminal_kmer_idx - 1);
    hmm_result.set3_by_event_kmer(terminal_event_idx, terminal_kmer_idx, score_d, score_u, score_l);

    return adaptive_banded_backtrack(hmm_result);
}

// Realign the read in signal space
void project_read(const faidx_t* fai,
                  const std::map<std::string, std::string>& fast5_path_map,
                  const Alphabet* alphabet,
                  const size_t k,
                  const std::vector<float>& levels,
                  const bam_hdr_t* hdr,
                  const bam1_t* record,
                  size_t read_idx,
                  int region_start,
                  int region_end)
{
    // Skip unmapped
    if((record->core.flag & BAM_FUNMAP) != 0) {
        return;
    }

    // Skip reads with hardclips
    const uint32_t *cigar = bam_get_cigar(record);
    if( bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP ||
        bam_cigar_op(cigar[record->core.n_cigar - 1]) == BAM_CHARD_CLIP ) 
     {
        return;
     }

    std::string read_name = bam_get_qname(record);

    // read move tag from bam
    uint8_t* mv_tag_ptr = bam_aux_get(record, "mv");
    if(mv_tag_ptr == NULL) {
        fprintf(stderr, "mv tag unavailable for %s\n", read_name.c_str());
        return;
    }

    uint32_t mv_len = bam_auxB_len(mv_tag_ptr);
    if(mv_len == 0) {
        fprintf(stderr, "mv tag is length 0 for %s\n", read_name.c_str());
        return;
    }

    uint32_t mv_stride = bam_auxB2i(mv_tag_ptr, 0);
    std::vector<uint32_t> mv;
    mv.reserve(mv_len - 1);
    for(int i = 1; i < mv_len; ++i) {
        mv.push_back(bam_auxB2i(mv_tag_ptr, i));
    }

    uint32_t mv_sum = 0;
    for(int i = 0; i < mv.size(); ++i) {
        mv_sum += mv[i];
    }
    
    // read trim parameter
    uint8_t* ts_tag_ptr = bam_aux_get(record, "ts");
    if(ts_tag_ptr == NULL) {
        fprintf(stderr, "ts tag unavailable for %s\n", read_name.c_str());
        return;
    }
    uint32_t mv_signal_start = bam_aux2i(ts_tag_ptr);

    // read shift parameters
    uint8_t* sm_tag_ptr = bam_aux_get(record, "sm");
    if(sm_tag_ptr == NULL) {
        fprintf(stderr, "sm tag unavailable for %s\n", read_name.c_str());
        return;
    }
    float basecaller_shift = bam_aux2f(sm_tag_ptr);

    // read scale parameter
    uint8_t* sd_tag_ptr = bam_aux_get(record, "sd");
    if(sd_tag_ptr == NULL) {
        fprintf(stderr, "sd tag unavailable for %s\n", read_name.c_str());
        return;
    }
    float basecaller_scale = bam_aux2f(sd_tag_ptr);

    uint8_t* f5_name_ptr = bam_aux_get(record, "f5");
    if(f5_name_ptr == NULL) {
        fprintf(stderr, "f5 tag unavailable for %s\n", read_name.c_str());
        return;
    }
    std::string f5_name(bam_aux2Z(f5_name_ptr));

    // If we've got here, everything in the bam is sufficient for projecting the signal

    // load the read sequence, on the original sequencing strand
    std::string read_sequence = bam_get_seq_str(record);
    bool is_rev = bam_is_rev(record);
    if(is_rev) {
        read_sequence = alphabet->reverse_complement(read_sequence);
    }

    // load the raw samples
    auto path_map_iter = fast5_path_map.find(f5_name);
    if(path_map_iter == fast5_path_map.end()) {
        fprintf(stderr, "could not find fast5 path for %s (expected filename %s not found in directory)\n", read_name.c_str(), f5_name.c_str());
        return;
    } 
    std::string fast5_file_path = path_map_iter->second;

    Fast5Data signal_data = Fast5Loader::load_read(fast5_file_path, read_name);

    // collect a vector with the starting sample for each base
    std::vector<int> start_samples;
    start_samples.reserve(read_sequence.length());

    // this is equivalent to what remora does, it takes the index of the non-zero elements
    // and multiplies by stride
    for(int i = 0; i < mv.size(); ++i) {
        if(mv[i] > 0) {
            start_samples.push_back(mv_signal_start + i * mv_stride);
        }
    }
    assert(start_samples.size() == read_sequence.length());

    // for convenience, to get the range for the last base
    start_samples.push_back(signal_data.rt.n);
    
    if(opt::verbose > 1) {
        fprintf(stderr, "[init] projecting %s (%zu samples, %zu bases)\n", read_name.c_str(), signal_data.rt.n, read_sequence.length());
    }

    /*
    fprintf(stderr, "[move parse] num moves: %zu sum of moves: %d read length: %zu ts: %d ss len: %zu basecaller scalings [%.1f %.1f]\n", 
        mv.size(), mv_sum, read_sequence.length(), mv_signal_start, start_samples.size(), 
        basecaller_shift, basecaller_scale);
    */

    // Extract the reference subsequence for the entire alignment
    int fetched_len = 0;
    int ref_offset = record->core.pos;
    std::string ref_name(hdr->target_name[record->core.tid]);
    std::string ref_seq = get_reference_region_ts(fai, ref_name.c_str(), ref_offset,
                                                  bam_endpos(record), &fetched_len);

    // convert to upper case
    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);

    // k from read pore model
    size_t strand_idx = 0;

    // If the reference sequence contains ambiguity codes
    // switch them to the lexicographically lowest base
    ref_seq = alphabet->disambiguate(ref_seq);
    std::string rc_ref_seq = alphabet->reverse_complement(ref_seq);

    // this was empirically determined by finding the value
    // that minimized the difference between scaled sample level
    // and the k-mer model
    int kmer_reg_offset = -5;
    
    // disallow trim and skips since we are aligning to the reads, which was the source
    // of initial move-based alignment
    float skip_pen = -INFINITY;
    float trim_pen = -INFINITY;

    // Refine signal mapping using DP.
    std::vector<float> expected_levels(read_sequence.length(), 0.0f);
    int64_t max_kmer_pos = read_sequence.length() - k;

    for(int64_t i = 0; i < read_sequence.length(); ++i) {
        int64_t kmer_pos = i + kmer_reg_offset;
        if(kmer_pos < 0) { kmer_pos = 0; }
        if(kmer_pos > max_kmer_pos) { kmer_pos = max_kmer_pos; }
        uint32_t rank = alphabet->kmer_rank(read_sequence.c_str() + kmer_pos, k);
        expected_levels[i] = levels[rank];
    }

    std::vector<float> scaled_samples(signal_data.rt.n, 0.0f);
    for(size_t i = 0; i < scaled_samples.size(); ++i) {
        double s = (signal_data.rt.raw[i] - basecaller_shift) / basecaller_scale;
        scaled_samples[i] = s;
    }

    std::vector<AlignedPair> refined_alignment = sample_dp(expected_levels, scaled_samples, start_samples, skip_pen, trim_pen);

    // reassign samples to bases
    for(int i = 0; i < start_samples.size() - 1; ++i) {
        start_samples[i] = -1;
    }

    int prev_idx = -1;
    for(int i = 0; i < refined_alignment.size(); ++i) {
        int curr_idx = refined_alignment[i].ref_pos;
        if(curr_idx != prev_idx) {
            start_samples[refined_alignment[i].ref_pos] = refined_alignment[i].read_pos;    
        }
        prev_idx = curr_idx;
    }
    
    float sum_abs_diff = 0.0;
    size_t samples_aligned = 0;
    size_t bases_aligned = 0;

    std::vector<AlignedSegment> aligned_segments = get_aligned_segments(record);
    for(size_t segment_idx = 0; segment_idx < aligned_segments.size(); ++segment_idx) {
        auto& segment = aligned_segments[segment_idx];
        
        // the read coordinates are wrt SEQ, which is always on the + strand
        // reverse these coordinates here so we can make sense of the signal data
        if(is_rev) {
            size_t l = read_sequence.length();
            for(size_t i = 0; i < segment.size(); ++i) {
                segment[i].read_pos = l - segment[i].read_pos - 1;
            }
        }

        for(size_t i = opt::alignment_end_trim; i < segment.size() - opt::alignment_end_trim; ++i) {
            auto& ap = segment[i];
            char ref_base = ref_seq[ap.ref_pos - ref_offset];
            char read_base = read_sequence[ap.read_pos];
            
            if(is_rev) {
                read_base = alphabet->complement(read_base);
            }

            // calculate signal summary
            float sample_sum = 0.0;
            int n_samples = 0;
            for(int j = start_samples[ap.read_pos]; j < start_samples[ap.read_pos+1]; ++j) {
                sample_sum += signal_data.rt.raw[j];
                n_samples += 1;
            }
        
            float mean = sample_sum / n_samples;
            float ont_norm_mean = (mean - basecaller_shift) / basecaller_scale;

            // get the reference kmer
            std::string kmer;
            kmer.reserve(k);
            if(!is_rev) {
                int kp = ap.ref_pos - ref_offset + kmer_reg_offset;
                if(kp < 0 || kp >= ref_seq.length()) {
                    continue;
                }
                kmer = ref_seq.substr(kp, k);
            } else {
                int seq_pos = ap.ref_pos - ref_offset;
                int rc_kp = rc_ref_seq.length() - seq_pos - 1 + kmer_reg_offset;
                if(rc_kp < 0 || rc_kp >= rc_ref_seq.length()) {
                    continue;
                }
                kmer = rc_ref_seq.substr(rc_kp, k);
            }

            if(kmer.length() != k) {
                continue;
            }
            uint32_t rank = alphabet->kmer_rank(kmer.c_str(), k);
            float model_level = levels[rank];
            float abs_diff = n_samples > 0 ? fabs(ont_norm_mean - model_level) : 0.0;
            sum_abs_diff += (n_samples * abs_diff);
            samples_aligned += n_samples;
            bases_aligned += start_samples[ap.read_pos] > -1;

            std::string samples_str = "-";
            std::stringstream samples_out;
            if(opt::write_samples && n_samples > 0) {
                samples_out.precision(4);
                for(int j = start_samples[ap.read_pos]; j < start_samples[ap.read_pos + 1]; ++j) {
                    double scaled_current = (signal_data.rt.raw[j] - basecaller_shift) / basecaller_scale;
                    samples_out << scaled_current << ",";
                }
                samples_str = samples_out.str();
            }

            #pragma omp critical
            printf("%s\t%d\t%c\t%s\t%c\t%d\t%c\t%.1f\t%.3f\t%s\t%.3f\t%d\t%d\t%s\t%.2f\n", 
                ref_name.c_str(), ap.ref_pos, ref_base, read_name.c_str(), !is_rev ? '+' : '-', ap.read_pos, read_base,
                mean, ont_norm_mean, kmer.c_str(), model_level, start_samples[ap.read_pos], start_samples[ap.read_pos + 1],
                samples_str.c_str(), sum_abs_diff);
        }
    }

    //printf("%s\t%d\t%.2f\t%.2f\t%.2lf\t%zu\t%zu\n", read_name.c_str(), kmer_reg_offset, skip_pen, trim_pen, sum_abs_diff, samples_aligned, bases_aligned);

    // cleanup 
    free(signal_data.rt.raw);
    signal_data.rt.raw = NULL;
}

void parse_project_signal_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'q': arg >> opt::min_mapping_quality; break;
            case 'm': arg >> opt::kmer_model_file; break;
            case 'd': arg >> opt::fast5_directory; break;
            case OPT_SAMPLES: opt::write_samples = true; break;
            case 'v': opt::verbose++; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << PROJECT_SIGNAL_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << PROJECT_SIGNAL_VERSION_MESSAGE;
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

    if(opt::genome_file.empty()) {
        std::cerr << SUBPROGRAM ": a --genome file must be provided\n";
        die = true;
    }

    if(opt::bam_file.empty()) {
        std::cerr << SUBPROGRAM ": a --bam file must be provided\n";
        die = true;
    }
    
    if(opt::kmer_model_file.empty()) {
        std::cerr << SUBPROGRAM ": a --kmer-model-file must be provided\n";
        die = true;
    }
    
    if(opt::fast5_directory.empty()) {
        std::cerr << SUBPROGRAM ": a --directory must be provided\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << PROJECT_SIGNAL_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

std::vector<float> load_ont_level_model(const std::string& filename, const Alphabet* alphabet, size_t& model_k)
{
    std::vector<float> levels;
    std::ifstream model_reader(filename);
    std::string model_line;
    while (getline(model_reader, model_line)) {
        if(model_line[0] == '#' || model_line.find("kmer") != -1) {
            fprintf(stderr, "Unexpected line in model file:\n");
            fprintf(stderr, "%s", model_line.c_str());
            fprintf(stderr, "Please use a model from ONT's repository: https://github.com/nanoporetech/kmer_models/");
            exit(EXIT_FAILURE);
        }

        std::stringstream parser(model_line);
        std::string kmer;
        float level;
        parser >> kmer >> level;

        model_k = kmer.length();
        size_t rank = alphabet->kmer_rank(kmer.c_str(), model_k);
        assert(rank == levels.size());
        levels.push_back(level);
    }
    return levels;
}

void find_fast5s(const std::string& path, std::map<std::string, std::string>& filenames)
{
    if (!is_directory(path)) {
        return;
    }

    auto dir_list = list_directory(path);
    for (const auto& fn : dir_list) {
        if(fn == "." or fn == "..") {
            continue;
        }

        std::string full_fn = path + "/" + fn;
        bool is_fast5 = full_fn.find(".fast5") != -1;
        if(is_fast5) { 
            // check if file exists in map
            if(filenames.find(fn) != filenames.end()) {
                std::string map_path = filenames[fn];
                fprintf(stderr, "error: %s already exists in directory tree:\n", fn.c_str());
                fprintf(stderr, "\tfirst path: %s\n", map_path.c_str());
                fprintf(stderr, "\tnew path: %s\n", full_fn.c_str());
                fprintf(stderr, "all fast5 filenames must be unique\n");
            }

            filenames[fn] = full_fn;
        } else {
            find_fast5s(full_fn, filenames);
        }
    }
}

int project_signal_main(int argc, char** argv)
{
    parse_project_signal_options(argc, argv);
    omp_set_num_threads(opt::num_threads);
    
    size_t model_k = 0;
    const Alphabet* alphabet = &gDNAAlphabet;

    // Load ONT k-mer model and pack it into a vector based on ranks over the alphabet
    std::vector<float> kmer_levels = load_ont_level_model(opt::kmer_model_file, alphabet, model_k);
    fprintf(stderr, "[init] Loaded %zu k-mer levels (k: %zu)\n", kmer_levels.size(), model_k);

    // Load paths to fast5 files
    std::map<std::string, std::string> fast5_path_map;
    find_fast5s(opt::fast5_directory, fast5_path_map);
    fprintf(stderr, "[init] Found %zu fast5 files\n", fast5_path_map.size());

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

#ifndef H5_HAVE_THREADSAFE
    if(opt::num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif
    printf("ref_name\tref_position\tref_base\tread_name\tstrand\tread_position\tread_base\traw_sample_mean\tbasecaller_scaled_mean\tref_kmer\tref_kmer_level\tsample_start\tsample_end\tsamples\tdebug\n");

    // the BamProcessor framework calls the input function with the
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    auto f = std::bind(project_read, std::ref(fai), std::ref(fast5_path_map), alphabet, model_k, std::ref(kmer_levels), _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.set_min_mapping_quality(opt::min_mapping_quality);

    // run
    processor.parallel_run(f);

    fai_destroy(fai);
    return EXIT_SUCCESS;
}
