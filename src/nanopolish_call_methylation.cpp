//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_call_methylation -- identify methylated bases
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <unordered_map>
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
#include <unistd.h>
#include <zlib.h>
#include "htslib/faidx.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_bam_processor.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_read_db.h"
#include "nanopolish_fast5_processor.h"
#include "fs_support.hpp"
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"
#include "minimap.h"
#include "bseq.h"
#include "mmpriv.h"

using namespace std::placeholders;

//
// Structs
//
struct OutputHandles
{
    FILE* site_writer = NULL;
    htsFile* bam_writer = NULL;

    void write_site_header() {
        // Write header
        fprintf(site_writer, "chromosome\tstrand\tstart\tend\tread_name\t"
                             "log_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\t"
                             "num_calling_strands\tnum_motifs\tsequence\n");
    }

    void close() {
        if(site_writer != stdout) {
            fclose(site_writer);
            site_writer = NULL;
        }

        if(bam_writer != NULL) {
            hts_close(bam_writer);
            bam_writer = NULL;
        }
    }
};

struct ScoredSite
{
    ScoredSite()
    {
        ll_unmethylated[0] = 0;
        ll_unmethylated[1] = 0;
        ll_methylated[0] = 0;
        ll_methylated[1] = 0;
        strands_scored = 0;
    }

    std::string chromosome;
    int start_position;
    int end_position;
    int n_motif;
    std::string sequence;

    // scores per strand
    double ll_unmethylated[2];
    double ll_methylated[2];
    int strands_scored;

    //
    static bool sort_by_position(const ScoredSite& a, const ScoredSite& b) { return a.start_position < b.start_position; }
};

struct WatchStatus
{
    WatchStatus() : timer("") {}
    ~WatchStatus() { fprintf(stderr, "\n"); } // end status line
    int processed_fast5s = 0;
    int total_fast5s = 0;
    Progress timer;

    void update(const std::string& message)
    {
        double elapsed_seconds = timer.get_elapsed_seconds();
        double seconds_per_file = processed_fast5s > 0 ? elapsed_seconds / processed_fast5s : 0.0;

        fprintf(stderr, "\33[2K\r[call-methylation] fast5 progress: %d/%d %.2lfs/fast5 [%s]", processed_fast5s, seconds_per_file, total_fast5s, message.c_str());
    }
};

//
const Alphabet* mtest_alphabet;

//
// Getopt
//
#define SUBPROGRAM "call-methylation"

static const char *CALL_METHYLATION_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *CALL_METHYLATION_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa --methylation cpg\n"
"Classify nucleotides as methylated or not.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -r, --reads=FILE                     the ONT reads are in fasta/fastq FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are calling methylation for is in fasta FILE\n"
"  -q, --methylation=STRING             the type of methylation (cpg,gpc,dam,dcm)\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --watch=DIR                      watch the sequencing run directory DIR and call methylation as data is generated\n"
"      --watch-write-bam                in watch mode, write the alignments for each fastq\n"
"  -c, --watch-process-total=TOTAL      in watch mode, there are TOTAL calling processes running against this directory\n"
"  -i, --watch-process-index=IDX        in watch mode, the index of this process is IDX\n"
"                                       the previous two options allow you to run multiple independent methylation\n"
"                                       calling processes against a single directory. Each process will only call\n"
"                                       files when X mod TOTAL == IDX, where X is the suffix of the fast5 file.\n"
"      --progress                       print out a progress message\n"
"  -K  --batchsize=NUM                  the batch size (default: 512)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string methylation_type = "cpg";
    static std::string models_fofn;
    static std::string region;
    static std::string motif_methylation_model_type = "reftrained";
    static std::string watch_dir;
    static int watch_write_bam;
    static int watch_process_total = 1;
    static int watch_process_index = 0;
    static int progress = 0;
    static int num_threads = 1;
    static int batch_size = 512;
    static int min_separation = 10;
    static int min_flank = 10;
    static int min_mapping_quality = 20;
}

static const char* shortopts = "r:b:g:t:w:m:K:q:c:i:vn";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PROGRESS, OPT_MIN_SEPARATION, OPT_WATCH_DIR, OPT_WATCH_WRITE_BAM };

static const struct option longopts[] = {
    { "verbose",              no_argument,       NULL, 'v' },
    { "reads",                required_argument, NULL, 'r' },
    { "bam",                  required_argument, NULL, 'b' },
    { "genome",               required_argument, NULL, 'g' },
    { "methylation",          required_argument, NULL, 'q' },
    { "window",               required_argument, NULL, 'w' },
    { "threads",              required_argument, NULL, 't' },
    { "models-fofn",          required_argument, NULL, 'm' },
    { "watch-process-total",  required_argument, NULL, 'c' },
    { "watch-process-index",  required_argument, NULL, 'i' },
    { "min-separation",       required_argument, NULL, OPT_MIN_SEPARATION },
    { "watch",                required_argument, NULL, OPT_WATCH_DIR },
    { "watch-write-bam",      no_argument,       NULL, OPT_WATCH_WRITE_BAM },
    { "progress",             no_argument,       NULL, OPT_PROGRESS },
    { "help",                 no_argument,       NULL, OPT_HELP },
    { "version",              no_argument,       NULL, OPT_VERSION },
    { "batchsize",            no_argument,       NULL, 'K' },
    { NULL, 0, NULL, 0 }
};

// Test motif sites in this read for methylation
void calculate_methylation_for_read(const OutputHandles& handles,
                                    SquiggleRead& sr,
                                    const faidx_t* fai,
                                    const bam_hdr_t* hdr,
                                    const bam1_t* record,
                                    size_t read_idx,
                                    int region_start,
                                    int region_end)
{
    std::string read_orientation = bam_is_rev(record) ? "-" : "+";

    // An output map from reference positions to scored motif sites
    std::map<int, ScoredSite> site_score_map;

    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {
        if(!sr.has_events_for_strand(strand_idx)) {
            continue;
        }

        size_t k = sr.get_model_k(strand_idx);

        // check if there is a motif model for this strand
        if(!PoreModelSet::has_model(sr.get_model_kit_name(strand_idx),
                                    opt::methylation_type,
                                    sr.get_model_strand_name(strand_idx),
                                    k))
        {
            continue;
        }

        // Build the event-to-reference map for this read from the bam record
        SequenceAlignmentRecord seq_align_record(record);
        EventAlignmentRecord event_align_record(&sr, strand_idx, seq_align_record);

        std::vector<double> site_scores;
        std::vector<int> site_starts;
        std::vector<int> site_ends;
        std::vector<int> site_count;


        std::string contig = hdr->target_name[record->core.tid];
        int ref_start_pos = record->core.pos;
        int ref_end_pos =  bam_endpos(record);

        // Extract the reference sequence for this region
        int fetched_len = 0;
        assert(ref_end_pos >= ref_start_pos);
        std::string ref_seq = get_reference_region_ts(fai, contig.c_str(), ref_start_pos,
                                                      ref_end_pos, &fetched_len);

        // Remove non-ACGT bases from this reference segment
        ref_seq = gDNAAlphabet.disambiguate(ref_seq);

        // Scan the sequence for motifs
        std::vector<int> motif_sites;
        assert(ref_seq.size() != 0);
        for(size_t i = 0; i < ref_seq.size() - 1; ++i) {
            if(mtest_alphabet->is_motif_match(ref_seq, i))
                motif_sites.push_back(i);
        }

        // Batch the motifs together into groups that are separated by some minimum distance
        std::vector<std::pair<int, int>> groups;

        size_t curr_idx = 0;
        while(curr_idx < motif_sites.size()) {
            // Find the endpoint of this group of sites
            size_t end_idx = curr_idx + 1;
            while(end_idx < motif_sites.size()) {
                if(motif_sites[end_idx] - motif_sites[end_idx - 1] > opt::min_separation)
                    break;
                end_idx += 1;
            }
            groups.push_back(std::make_pair(curr_idx, end_idx));
            curr_idx = end_idx;
        }

        for(size_t group_idx = 0; group_idx < groups.size(); ++group_idx) {

            size_t start_idx = groups[group_idx].first;
            size_t end_idx = groups[group_idx].second;

            // the coordinates on the reference substring for this group of sites
            int sub_start_pos = motif_sites[start_idx] - opt::min_flank;
            int sub_end_pos = motif_sites[end_idx - 1] + opt::min_flank;
            int span = motif_sites[end_idx - 1] - motif_sites[start_idx];

            // skip if too close to the start of the read alignment or
            // if the reference range is too large to efficiently call
            if(sub_start_pos <= opt::min_separation || span > 200) {
                continue;
            }

            std::string subseq = ref_seq.substr(sub_start_pos, sub_end_pos - sub_start_pos + 1);
            std::string rc_subseq = mtest_alphabet->reverse_complement(subseq);

            int calling_start = sub_start_pos + ref_start_pos;
            int calling_end = sub_end_pos + ref_start_pos;

            // using the reference-to-event map, look up the event indices for this segment
            int e1,e2;
            bool bounded = AlignmentDB::_find_by_ref_bounds(event_align_record.aligned_events,
                                                            calling_start,
                                                            calling_end,
                                                            e1,
                                                            e2);

            double ratio = fabs(e2 - e1) / (calling_start - calling_end);

            // Only process this region if the the read is aligned within the boundaries
            // and the span between the start/end is not unusually short
            if(!bounded || abs(e2 - e1) <= 10 || ratio > MAX_EVENT_TO_BP_RATIO) {
                continue;
            }

            uint32_t hmm_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;

            // Set up event data
            HMMInputData data;
            data.read = &sr;
            data.pore_model = sr.get_model(strand_idx, opt::methylation_type);
            data.strand = strand_idx;
            data.rc = event_align_record.rc;
            data.event_start_idx = e1;
            data.event_stop_idx = e2;
            data.event_stride = data.event_start_idx <= data.event_stop_idx ? 1 : -1;

            // Calculate the likelihood of the unmethylated sequence
            HMMInputSequence unmethylated(subseq, rc_subseq, mtest_alphabet);
            double unmethylated_score = profile_hmm_score(unmethylated, data, hmm_flags);

            // Methylate all motifs in the sequence and score again
            std::string m_subseq = mtest_alphabet->methylate(subseq);
            std::string rc_m_subseq = mtest_alphabet->reverse_complement(m_subseq);

            // Calculate the likelihood of the methylated sequence
            HMMInputSequence methylated(m_subseq, rc_m_subseq, mtest_alphabet);
            double methylated_score = profile_hmm_score(methylated, data, hmm_flags);

            // Aggregate score
            int start_position = motif_sites[start_idx] + ref_start_pos;
            auto iter = site_score_map.find(start_position);
            if(iter == site_score_map.end()) {
                // insert new score into the map
                ScoredSite ss;
                ss.chromosome = contig;
                ss.start_position = start_position;
                ss.end_position = motif_sites[end_idx - 1] + ref_start_pos;
                ss.n_motif = end_idx - start_idx;

                // extract the motif site(s) with a k-mers worth of surrounding context
                size_t site_output_start = motif_sites[start_idx] - k + 1;
                size_t site_output_end =  motif_sites[end_idx - 1] + k;
                ss.sequence = ref_seq.substr(site_output_start, site_output_end - site_output_start);

                // insert into the map
                iter = site_score_map.insert(std::make_pair(start_position, ss)).first;
            }

            // set strand-specific score
            // upon output below the strand scores will be summed
            iter->second.ll_unmethylated[strand_idx] = unmethylated_score;
            iter->second.ll_methylated[strand_idx] = methylated_score;
            iter->second.strands_scored += 1;
        } // for group
    } // for strands

    #pragma omp critical(call_methylation_write)
    {
        // write all sites for this read
        for(auto iter = site_score_map.begin(); iter != site_score_map.end(); ++iter) {

            const ScoredSite& ss = iter->second;
            double sum_ll_m = ss.ll_methylated[0] + ss.ll_methylated[1];
            double sum_ll_u = ss.ll_unmethylated[0] + ss.ll_unmethylated[1];
            double diff = sum_ll_m - sum_ll_u;

            // do not output if outside the window boundaries
            if((region_start != -1 && ss.start_position < region_start) ||
               (region_end != -1 && ss.end_position >= region_end)) {
                continue;
            }

            fprintf(handles.site_writer, "%s\t%s\t%d\t%d\t", ss.chromosome.c_str(), read_orientation.c_str(), ss.start_position, ss.end_position);
            fprintf(handles.site_writer, "%s\t%.2lf\t", sr.read_name.c_str(), diff);
            fprintf(handles.site_writer, "%.2lf\t%.2lf\t", sum_ll_m, sum_ll_u);
            fprintf(handles.site_writer, "%d\t%d\t%s\n", ss.strands_scored, ss.n_motif, ss.sequence.c_str());
        }
    }
}

void calculate_methylation_for_read_from_bam(const OutputHandles& handles,
                                             const ReadDB& read_db,
                                             const faidx_t* fai,
                                             const bam_hdr_t* hdr,
                                             const bam1_t* record,
                                             size_t read_idx,
                                             int region_start,
                                             int region_end)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);
    SquiggleRead sr(read_name, read_db);
    calculate_methylation_for_read(handles, sr, fai, hdr, record, read_idx, region_start, region_end);
}

void calculate_methylation_for_read_from_fast5(const OutputHandles& handles,
                                               const std::unordered_map<std::string, std::string>& sequence_map,
                                               const std::unordered_map<std::string, std::vector<bam1_t*>>& alignment_map,
                                               const faidx_t* fai,
                                               const bam_hdr_t* hdr,
                                               const Fast5Data& fast5_data)
{
    const std::string& read_name = fast5_data.read_name;
    const auto& a_iter = alignment_map.find(read_name);
    assert(a_iter != alignment_map.end());
    const std::vector<bam1_t*> bam_records = a_iter->second;

    if(bam_records.empty()) {
        // no alignment, skip
        return;
    }

    const auto& s_iter = sequence_map.find(read_name);
    if(s_iter == sequence_map.end()) {
        // no sequence, skip
        return;
    }

    //
    SquiggleRead sr(s_iter->second, fast5_data);

    //
    for(size_t i = 0; i < bam_records.size(); ++i) {
        calculate_methylation_for_read(handles, sr, fai, hdr, bam_records[i], -1, -1, -1);
    }
}

//
// Watch mode
//
struct FileBatch
{
    bool ready_to_call() { return !fast5_path.empty() && !fastq_path.empty() && !called; }

    std::string fast5_path;
    std::string fastq_path;
    bool called = false;
};

// init kseq reader
KSEQ_INIT(gzFile, gzread)

void bseq_destroy(mm_bseq1_t* s)
{
    free(s->name);
    s->name = NULL;

    free(s->seq);
    s->seq = NULL;

    if(s->qual) {
        free(s->qual);
        s->qual = NULL;
    }
    if(s->comment) {
        free(s->comment);
        s->comment = NULL;
    }
}

void populate_maps_from_fastq(OutputHandles& handles,
                              const std::string& fastq_filename,
                              std::unordered_map<std::string, std::string>& read_sequence_map,
                              std::unordered_map<std::string, std::vector<bam1_t*>>& read_alignment_map,
                              bam_hdr_t* hdr, mm_mapopt_t mopt, mm_idx_t* mi)
{
    // Read fastq file
    FILE* read_fp = fopen(fastq_filename.c_str(), "r");
    if(read_fp == NULL) {
        fprintf(stderr, "error: could not open %s for read\n", fastq_filename.c_str());
        exit(EXIT_FAILURE);
    }

    gzFile gz_read_fp = gzdopen(fileno(read_fp), "r");
    if(gz_read_fp == NULL) {
        fprintf(stderr, "error: could not open %s using gzdopen\n", fastq_filename.c_str());
        exit(EXIT_FAILURE);
    }

    int ret = 0;
    kseq_t* seq = kseq_init(gz_read_fp);
    while((ret = kseq_read(seq)) >= 0) {
        //
        read_sequence_map[seq->name.s] = seq->seq.s;
    }
    kseq_destroy(seq);
    seq = NULL;

    // get vector of names, so we can use parallel for
    std::vector<std::string> read_names;
    for(const auto& e : read_sequence_map) {
        read_names.push_back(e.first);

        // pre-populate alignment table
        read_alignment_map[e.first] = std::vector<bam1_t*>();
    }

    // initialize minimap2's thread storage
    std::vector<mm_tbuf_t*> tbuf(opt::num_threads, NULL);
    for(size_t i = 0; i < tbuf.size(); ++i) {
        tbuf[i] = mm_tbuf_init();
        assert(tbuf[i] != NULL);
    }

    #pragma omp parallel for
    for(size_t ri = 0; ri < read_names.size(); ++ri) {

        const std::string& read_name = read_names[ri];
        const std::string& sequence = read_sequence_map[read_name];

        //
        mm_reg1_t *reg;
        int j, i, n_reg;
        size_t thread_id = omp_get_thread_num();
        assert(thread_id < tbuf.size());

        reg = mm_map(mi, sequence.length(), sequence.c_str(), &n_reg, tbuf[thread_id], &mopt, 0); // get all hits for the query

        for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
            mm_reg1_t *r = &reg[j];
            assert(r->p); // with MM_F_CIGAR, this should not be NULL

            // build a sam record
            mm_bseq1_t bseq;
            bseq.name = strndup(read_name.c_str(), read_name.length());
            bseq.seq = strndup(sequence.c_str(), sequence.length());
            bseq.l_seq = sequence.length();
            bseq.qual = NULL;
            bseq.comment = NULL;

            kstring_t s = { 0, 0, NULL };
            mm_write_sam(&s, mi, &bseq, r, n_reg, reg);

            // convert record to bam
            bam1_t* record = bam_init1();
            int parse_ret = sam_parse1(&s, hdr, record);

            // if the write bam option is turned on, write the alignment to disk
            if(handles.bam_writer != NULL) {
                #pragma omp critical
                int write_ret = sam_write1(handles.bam_writer, hdr, record);
            }

            // store alignment
            if(record->core.qual >= opt::min_mapping_quality) {
                read_alignment_map[read_name].push_back(record);
            }
            bseq_destroy(&bseq);

            free(r->p);
            free(s.s);
            s = { 0, 0, NULL };
        }
        free(reg);
    }

    for(size_t i = 0; i < tbuf.size(); ++i) {
        mm_tbuf_destroy(tbuf[i]);
        tbuf[i] = NULL;
    }

    // clean up fastq reader
    gzclose(gz_read_fp);
    fclose(read_fp);
    seq = NULL;
}

//
bool process_batch(const std::string& basename, const FileBatch& batch, const faidx_t* fai, bam_hdr_t* hdr, mm_mapopt_t mopt, mm_idx_t* mi, WatchStatus& status)
{
    // Initialize writers
    OutputHandles handles;

    std::string calls_outname = basename + ".methylation_calls.tsv";
    handles.site_writer = fopen(calls_outname.c_str(), "w");
    handles.write_site_header();

    // optionall open a bam to write the alignments to
    if(opt::watch_write_bam) {
        std::string bam_outname = basename + ".bam";
        handles.bam_writer = hts_open(bam_outname.c_str(), "bw");
        int ret = sam_hdr_write(handles.bam_writer, hdr);
        if(ret != 0) {
            fprintf(stderr, "Error writing SAM header\n");
            exit(EXIT_FAILURE);
        }
    } else {
        handles.bam_writer = NULL;
    }

    // build indices:
    //  read_id -> read_sequence
    //  read_id -> primary alignment
    std::unordered_map<std::string, std::string> read_sequence_map;
    std::unordered_map<std::string, std::vector<bam1_t*>> read_alignment_map;
    status.update("aligning " + basename);
    populate_maps_from_fastq(handles, batch.fastq_path, read_sequence_map, read_alignment_map, hdr, mopt, mi);

    // use fast5 processor to iterate over the signal reads, and run our call-methylation function
    status.update("calling " + basename);
    Fast5Processor processor(batch.fast5_path, opt::num_threads);
    auto f = std::bind(calculate_methylation_for_read_from_fast5, std::ref(handles),
                                                                  std::ref(read_sequence_map),
                                                                  std::ref(read_alignment_map),
                                                                  fai,
                                                                  hdr,
                                                                  _1);

    // call.
    processor.parallel_run(f);

    // clean up bam records
    size_t count = 0;
    for(auto& x : read_alignment_map) {
        std::vector<bam1_t*> bam_records = x.second;
        for(size_t i = 0; i < bam_records.size(); ++i) {
            bam_destroy1(bam_records[i]);
            count += 1;
        }
        bam_records.clear();
    }
    handles.close();

    return true;
}

//
void run_watch_mode(const faidx_t* fai)
{
    WatchStatus status;

    // build sam/bam header structure by parsing the faidx,
    // converting it to a string then  string from fai file
    kstring_t header_str = { 0, 0, NULL };
    for (int i = 0; i < faidx_nseq(fai); ++i) {
        const char* tname = faidx_iseq(fai, i);
	    ksprintf(&header_str, "@SQ\tSN:%s\tLN:%d\n", tname, faidx_seq_len(fai, tname));
    }

    // Parse string to get header
    bam_hdr_t* hdr = sam_hdr_parse(header_str.l, header_str.s);
    free(header_str.s);
    header_str = { 0, 0, NULL };

    // load minimap2 index structures
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR; // perform alignment
    iopt.flag = 0; iopt.k = 15; // map-ont parameters

    // open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(opt::genome_file.c_str(), &iopt, 0);

    // read the first part of the index
    mm_idx_t *mi = mm_idx_reader_read(r, opt::num_threads);
    mm_mapopt_update(&mopt, mi);

    //  hacky way to reject multi-part indices, which I don't want to support now
    mm_idx_t *mi_done = mm_idx_reader_read(r, opt::num_threads);
    if(mi == NULL || mi_done != NULL) {
        fprintf(stderr, "Could not read minimap2 index\n");
        exit(EXIT_FAILURE);
    }

    // setup paths to watch
    std::string base_dir = opt::watch_dir;
    std::string fast5_dir = base_dir + "/fast5_pass";
    std::string fastq_dir = base_dir + "/fastq_pass";

    // map from basename -> FileBatch to track state of each file
    std::unordered_map<std::string, FileBatch> batches;

    // loop until we haven't seen new files for about 30 minutes
    const int MAX_ITERATIONS_UNTIL_STOP = 30;
    int iteration_count = 0;
    while(iteration_count < MAX_ITERATIONS_UNTIL_STOP) {
        status.update("checking for new files");

        // update file db with fast5s
        std::vector<std::string> fast5_files = list_directory(fast5_dir);

        for(const auto& fn : fast5_files) {
            if(ends_with(fn, ".fast5")) {
                std::string basename = strip_extension(fn, ".fast5");
                batches[basename].fast5_path = fast5_dir + "/" + fn;
            }
        }

        // update file db with fastqs
        std::vector<std::string> fastq_files = list_directory(fastq_dir);
        for(const auto& fn : fastq_files) {
            if(ends_with(fn, ".fastq")) {
                std::string basename = strip_extension(fn, ".fastq");
                batches[basename].fastq_path = fastq_dir + "/" + fn;
            }
        }

        // filter batches to the subset we want to handle in this process
        std::unordered_map<std::string, FileBatch>::iterator iter = batches.begin();
        while(iter != batches.end()) {
            size_t suffix_pos = iter->first.rfind("_");
            if(suffix_pos == std::string::npos) {
                fprintf(stderr, "Error: invalid suffix to file %s\n", iter->first.c_str());
                exit(EXIT_FAILURE);
            }

            std::stringstream parser(iter->first.substr(suffix_pos + 1));
            size_t file_numeric_suffix;
            parser >> file_numeric_suffix;
            if(file_numeric_suffix % opt::watch_process_total != opt::watch_process_index) {
                // remove file from batch
                iter = batches.erase(iter);
            } else {
                // keep file
                iter++;
            }
        }
        status.total_fast5s = batches.size();

        // iterate over collection and see which files need to be processed
        for(auto& e : batches) {
            if(e.second.ready_to_call()) {
                bool success = process_batch(e.first, e.second, fai, hdr, mopt, mi, status);
                e.second.called = success;
                status.processed_fast5s += 1;
                iteration_count = 0;
            }
        }

        status.update("Waiting for next batch - sleeping");
        iteration_count += 1;
        sleep(60);
    }

    // cleanup
    mm_idx_destroy(mi);
    mm_idx_reader_close(r);
    bam_hdr_destroy(hdr);
}

void run_from_bam(const faidx_t* fai)
{
    ReadDB read_db;
    read_db.load(opt::reads_file);

    // Initialize writers
    OutputHandles handles;
    handles.site_writer = stdout;
    handles.write_site_header();

    // the BamProcessor framework calls the input function with the
    // bam record, read index, etc passed as parameters
    // bind the other parameters the worker function needs here
    auto f = std::bind(calculate_methylation_for_read_from_bam, std::ref(handles), std::ref(read_db), fai, _1, _2, _3, _4, _5);
    BamProcessor processor(opt::bam_file, opt::region, opt::num_threads, opt::batch_size);
    processor.set_min_mapping_quality(opt::min_mapping_quality);
    processor.parallel_run(f);
}

void parse_call_methylation_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 'q': arg >> opt::methylation_type; break;
            case 'b': arg >> opt::bam_file; break;
            case '?': die = true; break;
            case 't': arg >> opt::num_threads; break;
            case 'm': arg >> opt::models_fofn; break;
            case 'w': arg >> opt::region; break;
            case 'v': opt::verbose++; break;
            case 'K': arg >> opt::batch_size; break;
            case 'c': arg >> opt::watch_process_total; break;
            case 'i': arg >> opt::watch_process_index; break;
            case OPT_MIN_SEPARATION: arg >> opt::min_separation; break;
            case OPT_WATCH_DIR: arg >> opt::watch_dir; break;
            case OPT_WATCH_WRITE_BAM: opt::watch_write_bam = true; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << CALL_METHYLATION_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CALL_METHYLATION_VERSION_MESSAGE;
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

    if(opt::watch_dir.empty()) {
        if(opt::reads_file.empty()) {
            std::cerr << SUBPROGRAM ": a --reads file must be provided\n";
            die = true;
        }

        if(opt::bam_file.empty()) {
            std::cerr << SUBPROGRAM ": a --bam file must be provided\n";
            die = true;
        }
    }

    if(opt::methylation_type.empty()) {
        std::cerr << SUBPROGRAM ": a --methylation type must be provided\n";
        die = true;
    } else {
        mtest_alphabet = get_alphabet_by_name(opt::methylation_type);
    }

    if(opt::watch_process_index >= opt::watch_process_total) {
        std::cerr << SUBPROGRAM ": invalid --watch-process-index (must be less than --watch-process-count)\n";
        die = true;
    }

    if(!opt::models_fofn.empty()) {
        // initialize the model set from the fofn
        PoreModelSet::initialize(opt::models_fofn);
    }

    if (die)
    {
        std::cout << "\n" << CALL_METHYLATION_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int call_methylation_main(int argc, char** argv)
{
    parse_call_methylation_options(argc, argv);

    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());
    if(fai == NULL) {
        fprintf(stderr, "Error: could not open genome file: %s\n", opt::genome_file.c_str());
        fprintf(stderr, "Please check the path is correct\n");
        exit(1);
    }

#ifndef H5_HAVE_THREADSAFE
    if(opt::num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif

    omp_set_num_threads(opt::num_threads);
    if(!opt::watch_dir.empty()) {
        run_watch_mode(fai);
    } else {
        run_from_bam(fai);
    }

    fai_destroy(fai);

    return EXIT_SUCCESS;
}

