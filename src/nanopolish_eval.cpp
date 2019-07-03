//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_eval -- framework for testing new
// models and methods
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <cmath>
#include <sys/time.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include <map>
#include <omp.h>
#include <getopt.h>
#include <cstddef>
#include <iterator>
#include <ostream>
#include "htslib/faidx.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_variant.h"
#include "nanopolish_haplotype.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_bam_processor.h"
#include "nanopolish_bam_utils.h"
#include "nanopolish_index.h"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"
#include "logger.hpp"

using namespace std::placeholders;

//
// structs
//

//
// Getopt
//
#define SUBPROGRAM "eval"

static const char *EVAL_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2019 Ontario Institute for Cancer Research\n";

static const char *EVAL_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa variants.vcf\n"
"Evaluate the performance of probabilistic models\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -r, --reads=FILE                     the reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome are in bam FILE\n"
"  -g, --genome=FILE                    the reference genome is in FILE\n"
"  -w, --window=STR                     only evaluate reads in the window STR (format: ctg:start-end)\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"  -a, --analysis-type=STR              type of analysis to perform (sequence-model or read-level)\n"
"  -l, --label=STR                      apply a user-defined label to the output\n"
"  -s, --substitions                    generate substitution errors\n"
"  -i, --insertions                     generate insertion errors\n"
"  -d, --deletions                      generate deletion errors\n"
"  -h, --homopolymer-indels             generate random indels in homopolymers\n"
"      --min-homopolymer-length=NUM     minimum length of homopolymer to consider (default: 5)\n"
"      --progress                       print out a progress message\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string region;

    static std::string label = "default";
    static std::string analysis_type = "sequence-model";
    static unsigned progress = 0;
    static unsigned num_threads = 1;
    static unsigned batch_size = 128;
    static size_t max_reads = -1;
    static int substitutions = 0;
    static int insertions = 0;
    static int deletions = 0;
    static int homopolymer_indels = 0;
    static int min_homopolymer_length = 5;
    static int min_flanking_sequence = 30;
}

static const char* shortopts = "r:b:g:t:w:m:e:a:vsdih";

enum { OPT_HELP = 1,
       OPT_VERSION,
       OPT_PROGRESS,
       OPT_MIN_HOMOPOLYMER_LENGTH,
       OPT_LOG_LEVEL
     };

static const struct option longopts[] = {
    { "verbose",                no_argument,       NULL, 'v' },
    { "reads",                  required_argument, NULL, 'r' },
    { "bam",                    required_argument, NULL, 'b' },
    { "genome",                 required_argument, NULL, 'g' },
    { "threads",                required_argument, NULL, 't' },
    { "window",                 required_argument, NULL, 'w' },
    { "max-reads",              required_argument, NULL, 'm' },
    { "analysis-type",          required_argument, NULL, 'a' },
    { "label",                  required_argument, NULL, 'l' },
    { "min-homopolymer-length", required_argument, NULL, OPT_MIN_HOMOPOLYMER_LENGTH},
    { "substitutions",          no_argument,       NULL, 's' },
    { "insertions",             no_argument,       NULL, 'd' },
    { "deletions",              no_argument,       NULL, 'i' },
    { "homopolymer-indels",     no_argument,       NULL, 'h' },
    { "progress",               no_argument,       NULL, OPT_PROGRESS },
    { "help",                   no_argument,       NULL, OPT_HELP },
    { "version",                no_argument,       NULL, OPT_VERSION },
    { "log-level",              required_argument, NULL, OPT_LOG_LEVEL },
    { NULL, 0, NULL, 0 }
};

void parse_eval_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'b': arg >> opt::bam_file; break;
            case 'w': arg >> opt::region; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'm': arg >> opt::max_reads; break;
            case 'a': arg >> opt::analysis_type; break;
            case 'l': arg >> opt::label; break;
            case 's': opt::substitutions = 1; break;
            case 'i': opt::insertions = 1; break;
            case 'd': opt::deletions = 1; break;
            case 'h': opt::homopolymer_indels = 1; break;
            case 'v': opt::verbose++; break;
            case OPT_MIN_HOMOPOLYMER_LENGTH: arg >> opt::min_homopolymer_length; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << EVAL_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << EVAL_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_LOG_LEVEL:
                logger::Logger::set_level_from_option(arg.str());
                break;
        }
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

    if (die) {
        std::cout << "\n" << EVAL_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

struct AccuracyStats
{
    int total_all = 0;
    int correct_all = 0;

    int total_indel = 0;
    int correct_indel = 0;

    int total_sub = 0;
    int correct_sub = 0;

    float get_accuracy_all() { return total_all > 0 ? (float)correct_all / total_all : 0.0f; }
    float get_accuracy_sub() { return total_sub > 0 ? (float)correct_sub / total_sub : 0.0f; }
    float get_accuracy_indel() { return total_indel > 0 ? (float)correct_indel / total_indel : 0.0f; }
};

AccuracyStats calculate_accuracy_stats(const std::vector<Variant> scored_variants)
{
    AccuracyStats stats;
    for(const Variant& v : scored_variants) {

        bool correct = v.quality > 0;

        stats.total_all += 1;
        stats.correct_all += correct;

        if(v.is_snp()) {
            stats.total_sub += 1;
            stats.correct_sub += correct;
        } else {
            stats.total_indel += 1;
            stats.correct_indel += correct;
        }
    }
    return stats;
}

//
//
//
class SequenceModel
{
    public:
        virtual float score(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags) const = 0;
        virtual std::string get_name(const PoreModel* pore_model) const = 0;
};

//
class ProfileHMMModel : public SequenceModel
{
    float score(const HMMInputSequence& sequence, const HMMInputData& data, const uint32_t flags) const
    {
        return profile_hmm_score(sequence, data, flags);
    }

    std::string get_name(const PoreModel* pore_model) const
    {
        return "np_hmm_" + pore_model->metadata.get_kit_name() + "_" + (char)(pore_model->k + 48) + "mer";
    }
};


//
//
//
class AnalysisType
{
    public:
        virtual void write_header() const = 0;
        virtual void process_scored_variants(const SquiggleRead& sr, const SequenceModel* model, const std::vector<Variant>& scored_variants) = 0;
        virtual void finalize() = 0;
};

class ReadVariantsAnalysis : public AnalysisType
{
    void write_header() const
    {
        fprintf(stdout, "ref_name\tref_position\tref_seq\talt_seq\tread_name\tmodel_name\tlog_likelihood_ratio\n");
    }

    void process_scored_variants(const SquiggleRead& sr, const SequenceModel* model, const std::vector<Variant>& scored_variants)
    {
        for(const Variant& v : scored_variants) {
           #pragma omp critical
            fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t%s\t%.2lf\n",
                    v.ref_name.c_str(), v.ref_position, v.ref_seq.c_str(), v.alt_seq.c_str(), sr.read_name.c_str(), model->get_name(sr.get_base_model(0)).c_str(), v.quality);
        }
    }

    void finalize() {} // intentionally does nothing
};

class ConsensusVariantsAnalysis : public AnalysisType
{
    public:
        ConsensusVariantsAnalysis() {}
        ConsensusVariantsAnalysis(const std::vector<Variant>& variants)
        {
            initialize(variants);
        }

        void initialize(const std::vector<Variant>& variants)
        {
            for(const auto& v : variants) {
                std::string key = v.key();
                variant_map[key] = v;
                variant_map[key].quality = 0.0f;
                depth_map[key] = 0;
            }
        }

        void write_header() const
        {
            fprintf(stdout, "ref_name\tref_position\tref_seq\talt_seq\tnum_reads\tlog_likelihood_ratio\n");
        }

        void process_scored_variants(const SquiggleRead& sr, const SequenceModel* model, const std::vector<Variant>& scored_variants)
        {
            for(const Variant& v : scored_variants) {
                std::string key = v.key();
                #pragma omp critical
                {
                    variant_map[key].quality += v.quality;
                    depth_map[key] += 1;
                }
            }
        }

        virtual void finalize()
        {
            for(const auto& kv : variant_map) {
                int depth = depth_map[kv.first];
                const Variant& v = kv.second;
                fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t%d\t%.2lf\n",
                        v.ref_name.c_str(), v.ref_position, v.ref_seq.c_str(), v.alt_seq.c_str(), opt::label.c_str(), depth, v.quality);
            }
        }

        std::map<std::string, Variant> variant_map;
        std::map<std::string, int> depth_map;
};

class ConsensusAccuracyAnalysis : public ConsensusVariantsAnalysis
{
    public:
        ConsensusAccuracyAnalysis(const std::vector<Variant>& variants)
        {
            initialize(variants);
        }

        void write_header() const
        {
            fprintf(stdout, "label\tcount_all\tcount_subs\tcount_indels\taccuracy_all\taccuracy_subs\taccuracy_indels\n");
        }

        void finalize()
        {
            std::vector<Variant> variants;
            for(const auto& kv : variant_map) {
                variants.push_back(kv.second);
            }

            AccuracyStats accuracy_stats = calculate_accuracy_stats(variants);

            #pragma omp critical
            fprintf(stdout, "%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", opt::label.c_str(),
                accuracy_stats.total_all, accuracy_stats.total_sub, accuracy_stats.total_indel,
                accuracy_stats.get_accuracy_all(), accuracy_stats.get_accuracy_sub(), accuracy_stats.get_accuracy_indel());
        }
};

class ReadAccuracyAnalysis : public AnalysisType
{
    void write_header() const
    {
        fprintf(stdout, "read_name\tlabel\tcurrent_range\tmedian_sd\tshift\tscale\tvar\tread_rate\tevents_per_base\t"
                        "extra_event_frequency\tskip_frequency\tstall_frequency\taccuracy_all\taccuracy_subs\taccuracy_indels\n");
    }

    void process_scored_variants(const SquiggleRead& sr, const SequenceModel* model, const std::vector<Variant>& scored_variants)
    {
        // read metadata
        int strand_idx = 0;

        // speed
        double read_duration = sr.events[strand_idx].back().start_time - sr.events[strand_idx].front().start_time;
        int read_length = sr.read_sequence.length();
        double read_rate = read_length / read_duration;

        RateMetrics rate_metrics = sr.calculate_rate_metrics(strand_idx);
        SNRMetrics snr_metrics = sr.calculate_snr_metrics(strand_idx);
        AccuracyStats accuracy_stats = calculate_accuracy_stats(scored_variants);

        if(accuracy_stats.total_all > 0) {
            #pragma omp critical
            fprintf(stdout, "%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.4f\n",
                sr.read_name.c_str(), opt::label.c_str(), snr_metrics.current_range, snr_metrics.median_sd,
                sr.scalings[strand_idx].shift, sr.scalings[strand_idx].scale, sr.scalings[strand_idx].var,
                read_rate, sr.events_per_base[strand_idx], rate_metrics.extra_event_frequency, rate_metrics.skip_frequency, rate_metrics.stall_frequency,
                accuracy_stats.get_accuracy_all(), accuracy_stats.get_accuracy_sub(), accuracy_stats.get_accuracy_indel());
        }
    }

    void finalize() {} // do nothing
};

void generate_random_substitutions(const Haplotype& reference, const int num_variants, std::vector<Variant>& variants)
{
    int start = reference.get_reference_position();
    int length = reference.get_sequence().length();
    int adjusted_length = length - 2 * opt::min_flanking_sequence;

    for(size_t i = 0; i < num_variants; ++i) {
        int variant_position = start + opt::min_flanking_sequence + rand() % adjusted_length;
        std::string refseq = reference.substr_by_reference(variant_position, variant_position).get_sequence();
        std::string possible;
        for(size_t bi = 0; bi < 4; bi++) {
            char b = "ACGT"[bi];
            if(b != refseq[0]) {
                possible.append(1, b);
            }
        }
        std::string altseq(1, possible[rand() % possible.length()]);

        Variant v;
        v.ref_name = reference.get_reference_name();
        v.ref_position = variant_position;
        v.ref_seq = refseq;
        v.alt_seq = altseq;
        variants.push_back(v);
    }
}

void generate_random_insertions(const Haplotype& reference, const int num_variants, std::vector<Variant>& variants)
{
    int start = reference.get_reference_position();
    int length = reference.get_sequence().length();
    int adjusted_length = length - 2 * opt::min_flanking_sequence;

    for(size_t i = 0; i < num_variants; ++i) {
        int variant_position = start + opt::min_flanking_sequence + rand() % adjusted_length;
        std::string refseq = reference.substr_by_reference(variant_position, variant_position).get_sequence();
        std::string possible = "ACGT";
        std::string altseq = refseq + possible[rand() % possible.length()];

        Variant v;
        v.ref_name = reference.get_reference_name();
        v.ref_position = variant_position;
        v.ref_seq = refseq;
        v.alt_seq = altseq;
        variants.push_back(v);
    }
}

void generate_random_deletions(const Haplotype& reference, const int num_variants, std::vector<Variant>& variants)
{
    int start = reference.get_reference_position();
    int length = reference.get_sequence().length();
    int adjusted_length = length - 2 * opt::min_flanking_sequence;

    for(size_t i = 0; i < num_variants; ++i) {
        int variant_position = start + opt::min_flanking_sequence + rand() % adjusted_length;
        std::string refseq = reference.substr_by_reference(variant_position, variant_position + 1).get_sequence();
        std::string altseq = refseq.substr(0, 1);

        Variant v;
        v.ref_name = reference.get_reference_name();
        v.ref_position = variant_position;
        v.ref_seq = refseq;
        v.alt_seq = altseq;
        variants.push_back(v);
    }
}

void generate_random_homopolymer_errors(const Haplotype& reference, const int min_hp_length, float p_deletion, const int num_variants, std::vector<Variant>& variants)
{
    int start = reference.get_reference_position();
    int length = reference.get_sequence().length();
    int adjusted_length = length - 2 * opt::min_flanking_sequence;
    const std::string& reference_sequence = reference.get_sequence();

    // scan for homopolymers
    std::vector<Variant> candidate_variants;
    size_t i = opt::min_flanking_sequence;
    size_t ref_end = reference_sequence.length() - opt::min_flanking_sequence;
    while(i < ref_end) {

        // start a new homopolymer
        char hp_base = reference_sequence[i];
        size_t hap_hp_start = i;

        // extend until a different base is found
        while(i < ref_end && reference_sequence[i] == hp_base) i++;

        size_t ref_hp_start = hap_hp_start + reference.get_reference_position();
        size_t hap_hp_end = i;
        size_t hp_length = hap_hp_end - hap_hp_start;

        // disallow variants in the flanking region
        if(i >= ref_end || hp_length < min_hp_length) {
            continue;
        }

        // Construct the variant
        std::string refseq = reference.substr_by_reference(ref_hp_start, ref_hp_start + hp_length - 1).get_sequence();
        assert(refseq.length() == hp_length);
        std::string altseq = refseq;
        if( (float)rand() / RAND_MAX < p_deletion) {
            // deletion
            altseq = altseq.substr(0, altseq.length() - 1);
        } else {
            // insertion
            altseq = altseq.append(1, hp_base);
        }
        Variant v;
        v.ref_name = reference.get_reference_name();
        v.ref_position = ref_hp_start;
        v.ref_seq = refseq;
        v.alt_seq = altseq;
        candidate_variants.push_back(v);
    }

    // shuffle and take the top N
    std::random_shuffle(candidate_variants.begin(), candidate_variants.end());
    candidate_variants.resize(num_variants);
    variants.insert(variants.end(), candidate_variants.begin(), candidate_variants.end());
    return;
}


std::vector<Variant> score_variants_single_read(SquiggleRead& sr,
                                                const faidx_t* fai,
                                                const std::vector<Variant>& variants,
                                                const bam_hdr_t* hdr,
                                                const bam1_t* record,
                                                const SequenceModel* model)
{
    //
    std::vector<Variant> out_variants;
    uint32_t alignment_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;
    size_t strand_idx = 0;

    std::string ref_name = hdr->target_name[record->core.tid];
    int alignment_start_pos = record->core.pos;
    int alignment_end_pos = bam_endpos(record);

    // skip if region too short
    if(alignment_end_pos - alignment_start_pos < 1000) {
        return out_variants;
    }

    // load reference sequence
    int reference_length;
    std::string reference_seq =
        get_reference_region_ts(fai, ref_name.c_str(), alignment_start_pos, alignment_end_pos, &reference_length);
    Haplotype reference_haplotype(ref_name, alignment_start_pos, reference_seq);

    SequenceAlignmentRecord seq_align_record(record);
    EventAlignmentRecord event_align_record(&sr, strand_idx, seq_align_record);


    for(const Variant& v : variants) {

        int calling_start = v.ref_position - opt::min_flanking_sequence;
        int calling_end = v.ref_position + opt::min_flanking_sequence;

        if(v.ref_position <= alignment_start_pos || v.ref_position >= alignment_end_pos) {
            continue;
        }

        HMMInputData data;
        data.read = event_align_record.sr;
        data.strand = event_align_record.strand;
        data.rc = event_align_record.rc;
        data.event_stride = event_align_record.stride;
        data.pore_model = sr.get_base_model(event_align_record.strand);
        int e1,e2;
        bool bounded = AlignmentDB::_find_by_ref_bounds(event_align_record.aligned_events,
                                                        calling_start,
                                                        calling_end,
                                                        e1,
                                                        e2);

        // The events of this read do not span the calling window, skip
        if(!bounded || fabs(e2 - e1) / (calling_start - calling_end) > MAX_EVENT_TO_BP_RATIO) {
            continue;
        }

        data.event_start_idx = e1;
        data.event_stop_idx = e2;

        Haplotype calling_haplotype =
            reference_haplotype.substr_by_reference(calling_start, calling_end);

        //fprintf(stderr, "ref: %s\n", calling_haplotype.get_sequence().c_str());
        double ref_score = model->score(calling_haplotype.get_sequence(), data, alignment_flags);
        bool good_haplotype = calling_haplotype.apply_variant(v);
        if(good_haplotype) {
            double alt_score = model->score(calling_haplotype.get_sequence(), data, alignment_flags);
            out_variants.push_back(v);
            out_variants.back().quality = ref_score - alt_score;
            /*
            #pragma omp critical
            fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t%s\t%.2lf\t%.2lf\t%.2lf\n",
                ref_name.c_str(), v.ref_position, v.ref_seq.c_str(), v.alt_seq.c_str(), data.read->read_name.c_str(), model->get_name(data).c_str(), ref_score, alt_score, ref_score - alt_score);
            */
        }
    }
    return out_variants;
}

//
void eval_single_read(const ReadDB& read_db,
                      const faidx_t* fai,
                      const std::vector<Variant>& variants,
                      AnalysisType* analysis_type,
                      const bam_hdr_t* hdr,
                      const bam1_t* record,
                      size_t read_idx,
                      int region_start,
                      int region_end)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);
    SquiggleRead sr(read_name, read_db);

    // skip if no events
    if(!sr.has_events_for_strand(0)) {
        return;
    }

    SequenceModel* model = new ProfileHMMModel;

    std::vector<Variant> scored_variants = score_variants_single_read(sr, fai, variants, hdr, record, model);
    std::string ref_name = hdr->target_name[record->core.tid];
    analysis_type->process_scored_variants(sr, model, scored_variants);

    delete model;
}

int eval_main(int argc, char** argv)
{
    parse_eval_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    std::string contig;
    int start_base;
    int end_base;
    if(!opt::region.empty()) {
        parse_region_string(opt::region, contig, start_base, end_base);
    } else {
        fprintf(stderr, "Error: a region string must be provided");
        exit(EXIT_FAILURE);
    }

    // Generate variants across the reference region
    std::vector<Variant> variants;

    // load reference sequence
    int reference_length;
    std::string reference_seq =
        get_reference_region_ts(fai, contig.c_str(), start_base, end_base, &reference_length);

    Haplotype reference_haplotype(contig, start_base, reference_seq);

    // Create random variants in this reference sequence
    float variant_rate = 0.01;
    int num_variants = variant_rate * reference_length;
    if(opt::substitutions)
        generate_random_substitutions(reference_haplotype, num_variants, variants);
    if(opt::deletions)
        generate_random_deletions(reference_haplotype, num_variants, variants);
    if(opt::insertions)
        generate_random_insertions(reference_haplotype, num_variants, variants);
    if(opt::homopolymer_indels)
        generate_random_homopolymer_errors(reference_haplotype, opt::min_homopolymer_length, 0.5, num_variants, variants);

    AnalysisType* analysis_type;
    if(opt::analysis_type == "read-variants") {
        analysis_type = new ReadVariantsAnalysis;
    } else if(opt::analysis_type == "read-accuracy") {
        analysis_type = new ReadAccuracyAnalysis;
    } else if(opt::analysis_type == "consensus-variants") {
        analysis_type = new ConsensusVariantsAnalysis(variants);
    } else if(opt::analysis_type == "consensus-accuracy") {
        analysis_type = new ConsensusAccuracyAnalysis(variants);
    } else {
        fprintf(stderr, "Error: unknown analysis type %s\n", opt::analysis_type.c_str());
        exit(EXIT_FAILURE);
    }

    analysis_type->write_header();

    // the BamProcessor framework calls the input function with the
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads);
    processor.set_max_reads(opt::max_reads);
    processor.set_min_mapping_quality(20);

    auto f = std::bind(eval_single_read, std::ref(read_db), std::ref(fai), std::ref(variants), analysis_type, _1, _2, _3, _4, _5);
    processor.parallel_run(f);

    analysis_type->finalize();

    fai_destroy(fai);
    return EXIT_SUCCESS;
}
