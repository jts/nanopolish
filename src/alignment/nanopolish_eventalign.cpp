//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_consensus.cpp -- entry point to consensus functions
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
#include <set>
#include <omp.h>
#include <getopt.h>
#include <iterator>
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
#include "H5pubconf.h"
#include "profiler.h"
#include "progress.h"

//
// Getopt
//
#define SUBPROGRAM "eventalign"

static const char *EVENTALIGN_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *EVENTALIGN_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Align nanopore events to reference k-mers\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"      --sam                            write output in SAM format\n"
"  -w, --window=STR                     compute the consensus for window STR (format: ctg:start_id-end_id)\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"      --scale-events                   scale events to the model, rather than vice-versa\n"
"      --progress                       print out a progress message\n"
"  -n, --print-read-names               print read names instead of indexes\n"
"      --summary=FILE                   summarize the alignment of each read/strand in FILE\n"
"      --samples                        write the raw samples for the event to the tsv output\n"
"      --models-fofn=FILE               read alternative k-mer models from FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string region;
    static std::string summary_file;
    static std::string models_fofn;
    static int output_sam = 0;
    static int progress = 0;
    static int num_threads = 1;
    static int scale_events = 0;
    static int batch_size = 128;
    static bool print_read_names;
    static bool full_output;
    static bool write_samples = false;
}

static const char* shortopts = "r:b:g:t:w:vn";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PROGRESS, OPT_SAM, OPT_SUMMARY, OPT_SCALE_EVENTS, OPT_MODELS_FOFN, OPT_SAMPLES };

static const struct option longopts[] = {
    { "verbose",          no_argument,       NULL, 'v' },
    { "reads",            required_argument, NULL, 'r' },
    { "bam",              required_argument, NULL, 'b' },
    { "genome",           required_argument, NULL, 'g' },
    { "window",           required_argument, NULL, 'w' },
    { "threads",          required_argument, NULL, 't' },
    { "summary",          required_argument, NULL, OPT_SUMMARY },
    { "models-fofn",      required_argument, NULL, OPT_MODELS_FOFN },
    { "print-read-names", no_argument,       NULL, 'n' },
    { "samples",          no_argument,       NULL, OPT_SAMPLES },
    { "scale-events",     no_argument,       NULL, OPT_SCALE_EVENTS },
    { "sam",              no_argument,       NULL, OPT_SAM },
    { "progress",         no_argument,       NULL, OPT_PROGRESS },
    { "help",             no_argument,       NULL, OPT_HELP },
    { "version",          no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

// convenience wrapper for the two output modes
struct EventalignWriter
{
    FILE* tsv_fp;
    htsFile* sam_fp;
    FILE* summary_fp;
};

// Summarize the event alignment for a read strand
struct EventalignSummary
{
    EventalignSummary() {
        num_events = 0;
        num_steps = 0;
        num_stays = 0;
        num_skips = 0;
        sum_z_score = 0;
        sum_duration = 0;
        alignment_edit_distance = 0;
        reference_span = 0;
    }

    int num_events;
    int num_steps;
    int num_stays;
    int num_skips;

    double sum_duration;
    double sum_z_score;
    int alignment_edit_distance;
    int reference_span;
};

//
const PoreModel* EventAlignmentParameters::get_model() const 
{
    if(this->alphabet == "") {
        return this->sr->get_base_model(this->strand_idx);
    } else {
        return this->sr->get_model(this->strand_idx, this->alphabet);
    }
}

// Modify the aligned_pairs vector to ensure the highest read position
// does not exceed max_kmer
void trim_aligned_pairs_to_kmer(std::vector<AlignedPair>& aligned_pairs, int max_kmer_idx)
{
    int idx = aligned_pairs.size() - 1;
    while(idx >= 0 && aligned_pairs[idx].read_pos > max_kmer_idx)
        idx -= 1;

    if(idx < 0)
        aligned_pairs.clear(); // no valid data
    else
        aligned_pairs.resize(idx + 1);
}

// Modify the aligned_pairs vector to ensure there are no alignments
// outside of the given reference coordinates
void trim_aligned_pairs_to_ref_region(std::vector<AlignedPair>& aligned_pairs, int ref_start, int ref_end)
{
    std::vector<AlignedPair> trimmed;
    for(size_t i = 0; i < aligned_pairs.size(); ++i) {
        if(aligned_pairs[i].ref_pos >= ref_start && 
           aligned_pairs[i].ref_pos <= ref_end) {
            trimmed.push_back(aligned_pairs[i]);
        }
    }
    
    aligned_pairs.swap(trimmed);
}

// Returns the index into the aligned_pairs vector that has the highest ref_pos
// that is not greater than ref_pos_max. It starts the search at pair_idx
int get_end_pair(const std::vector<AlignedPair>& aligned_pairs, int ref_pos_max, int pair_idx)
{
    while(pair_idx < (int)aligned_pairs.size()) {
        if(aligned_pairs[pair_idx].ref_pos > ref_pos_max)
            return pair_idx - 1;
        pair_idx += 1;
    }
    
    return aligned_pairs.size() - 1;
}

// get the specified reference region, threadsafe
std::string get_reference_region_ts(const faidx_t* fai, const char* ref_name, int start, int end, int* fetched_len)
{

    // faidx_fetch_seq is not threadsafe
    char* cref_seq;
    #pragma omp critical
    cref_seq = faidx_fetch_seq(fai, ref_name, start, end, fetched_len);
    
    assert(cref_seq != NULL);

    std::string out(cref_seq);
    free(cref_seq);
    return out;
}

//
//
//

void emit_tsv_header(FILE* fp)
{
    fprintf(fp, "%s\t%s\t%s\t%s\t%s\t", "contig", "position", "reference_kmer",
            (not opt::print_read_names? "read_index" : "read_name"), "strand");
    fprintf(fp, "%s\t%s\t%s\t%s\t", "event_index", "event_level_mean", "event_stdv", "event_length");
    fprintf(fp, "%s\t%s\t%s\t%s", "model_kmer", "model_mean", "model_stdv", "standardized_level");

    if(opt::write_samples) {
        fprintf(fp, "\t%s", "samples");
    }
    fprintf(fp, "\n");
}

void emit_sam_header(samFile* fp, const bam_hdr_t* hdr)
{
    sam_hdr_write(fp, hdr);
}

std::string cigar_ops_to_string(const std::vector<uint32_t>& ops)
{
    std::stringstream ss;
    for(size_t i = 0; i < ops.size(); ++i) {
        ss << bam_cigar_oplen(ops[i]);
        ss << BAM_CIGAR_STR[bam_cigar_op(ops[i])];
    }
    return ss.str();
}

std::vector<uint32_t> event_alignment_to_cigar(const std::vector<EventAlignment>& alignments)
{
    std::vector<uint32_t> out;

    // add a softclip tag to account for unaligned events at the beginning/end of the read
    if(alignments[0].event_idx > 0) {
        out.push_back(alignments[0].event_idx << BAM_CIGAR_SHIFT | BAM_CSOFT_CLIP);
    }

    // we always start with a match
    out.push_back(1 << BAM_CIGAR_SHIFT | BAM_CMATCH);

    int prev_r_idx = alignments[0].ref_position;
    int prev_e_idx = alignments[0].event_idx;
    size_t ai = 1;

    while(ai < alignments.size()) {

        int r_idx = alignments[ai].ref_position;
        int e_idx = alignments[ai].event_idx;

        int r_step = abs(r_idx - prev_r_idx);
        int e_step = abs(e_idx - prev_e_idx);

        uint32_t incoming;
        if(r_step == 1 && e_step == 1) {

            // regular match
            incoming = 1 << BAM_CIGAR_SHIFT;
            incoming |= BAM_CMATCH;

        } else if(r_step > 1) {
            assert(e_step == 1);
            // reference jump of more than 1, this is how deletions are represented
            // we push the deletion onto the output then start a new match
            incoming = (r_step - 1) << BAM_CIGAR_SHIFT;
            incoming |= BAM_CDEL;
            out.push_back(incoming);
            
            incoming = 1 << BAM_CIGAR_SHIFT;
            incoming |= BAM_CMATCH;
        } else {
            assert(e_step == 1 && r_step == 0);
            incoming = 1 << BAM_CIGAR_SHIFT;
            incoming |= BAM_CINS;
        }

        // If the operation matches the previous, extend the length
        // otherwise append a new op
        if(bam_cigar_op(out.back()) == bam_cigar_op(incoming)) {
            uint32_t sum = bam_cigar_oplen(out.back()) + 
                           bam_cigar_oplen(incoming);
            out.back() = sum << BAM_CIGAR_SHIFT | bam_cigar_op(incoming);
        } else {
            out.push_back(incoming);
        }

        prev_r_idx = r_idx;
        prev_e_idx = e_idx;
        ai++;
    }
    return out;
}

void emit_event_alignment_sam(htsFile* fp,
                              const SquiggleRead& sr,
                              const bam_hdr_t* base_hdr,
                              const bam1_t* base_record, 
                              const std::vector<EventAlignment>& alignments)
{
    if(alignments.empty())
        return;
    bam1_t* event_record = bam_init1();
    
    // Variable-length data
    std::string qname = sr.read_name + (alignments.front().strand_idx == 0 ? ".template" : ".complement");

    // basic stats
    event_record->core.tid = base_record->core.tid;
    event_record->core.pos = alignments.front().ref_position;
    event_record->core.qual = base_record->core.qual;
    event_record->core.l_qname = qname.length() + 1; // must be null-terminated

    event_record->core.flag = alignments.front().rc ? 16 : 0;

    event_record->core.l_qseq = 0;
    
    event_record->core.mtid = -1;
    event_record->core.mpos = -1;
    event_record->core.isize = 0;

    std::vector<uint32_t> cigar = event_alignment_to_cigar(alignments);
    event_record->core.n_cigar = cigar.size();

    // calculate length of incoming data
    event_record->m_data = event_record->core.l_qname + // query name
                           event_record->core.n_cigar * 4 + // 4 bytes per cigar op
                           event_record->core.l_qseq + // query seq
                           event_record->core.l_qseq; // query quality
        
    // nothing copied yet
    event_record->l_data = 0;
    
    // allocate data
    event_record->data = (uint8_t*)malloc(event_record->m_data);

    // copy q name
    assert(event_record->core.l_qname <= event_record->m_data);
    strncpy(bam_get_qname(event_record), 
            qname.c_str(),
            event_record->core.l_qname);
    event_record->l_data += event_record->core.l_qname;
    
    // cigar
    assert(event_record->l_data + event_record->core.n_cigar * 4 <= event_record->m_data);
    memcpy(bam_get_cigar(event_record), 
           &cigar[0],
           event_record->core.n_cigar * 4);
    event_record->l_data += event_record->core.n_cigar * 4;

    // no copy for seq and qual
    assert(event_record->l_data <= event_record->m_data);

    int stride = alignments.front().event_idx < alignments.back().event_idx ? 1 : -1;
    bam_aux_append(event_record, "ES", 'i', 4, reinterpret_cast<uint8_t*>(&stride));

    sam_write1(fp, base_hdr, event_record);
    bam_destroy1(event_record); // automatically frees malloc'd segment
}

void emit_event_alignment_tsv(FILE* fp,
                              const SquiggleRead& sr,
                              uint32_t strand_idx,
                              const EventAlignmentParameters& params,
                              const std::vector<EventAlignment>& alignments)
{
    assert(params.alphabet == "");
    const PoreModel* pore_model = params.get_model();
    uint32_t k = pore_model->k;
    for(size_t i = 0; i < alignments.size(); ++i) {

        const EventAlignment& ea = alignments[i];

        // basic information
        if (not opt::print_read_names)
        {
            fprintf(fp, "%s\t%d\t%s\t%zu\t%c\t",
                    ea.ref_name.c_str(),
                    ea.ref_position,
                    ea.ref_kmer.c_str(),
                    ea.read_idx,
                    "tc"[ea.strand_idx]);
        }
        else
        {
            fprintf(fp, "%s\t%d\t%s\t%s\t%c\t",
                    ea.ref_name.c_str(),
                    ea.ref_position,
                    ea.ref_kmer.c_str(),
                    sr.read_name.c_str(),
                    "tc"[ea.strand_idx]);
        }

        // event information
        float event_mean = sr.get_unscaled_level(ea.event_idx, ea.strand_idx);
        float event_stdv = sr.get_stdv(ea.event_idx, ea.strand_idx);
        float event_duration = sr.get_duration(ea.event_idx, ea.strand_idx);
        uint32_t rank = pore_model->pmalphabet->kmer_rank(ea.model_kmer.c_str(), k);
        float model_mean = 0.0;
        float model_stdv = 0.0;

        if(opt::scale_events) {

            // scale reads to the model
            event_mean = sr.get_fully_scaled_level(ea.event_idx, ea.strand_idx);

            // unscaled model parameters
            if(ea.hmm_state != 'B') {
                PoreModelStateParams model = pore_model->get_parameters(rank);
                model_mean = model.level_mean;
                model_stdv = model.level_stdv;
            }
        } else {

            // scale model to the reads
            if(ea.hmm_state != 'B') {
                GaussianParameters model = sr.get_scaled_gaussian_from_pore_model_state(*pore_model, ea.strand_idx, rank);
                model_mean = model.mean;
                model_stdv = model.stdv;
            }
        }

        float standard_level = (event_mean - model_mean) / (sqrt(sr.scalings[ea.strand_idx].var) * model_stdv);
        fprintf(fp, "%d\t%.2lf\t%.3lf\t%.5lf\t", ea.event_idx, event_mean, event_stdv, event_duration);
        fprintf(fp, "%s\t%.2lf\t%.2lf\t%.2lf", ea.model_kmer.c_str(),
                                               model_mean,
                                               model_stdv,
                                               standard_level);

        if(opt::write_samples) {
            std::vector<float> samples = sr.get_scaled_samples_for_event(ea.strand_idx, ea.event_idx);
            std::stringstream sample_ss;
            std::copy(samples.begin(), samples.end(), std::ostream_iterator<float>(sample_ss, ","));

            // remove training comma
            std::string sample_str = sample_ss.str();
            sample_str.resize(sample_str.size() - 1);
            fprintf(fp, "\t%s", sample_str.c_str());
        }
        fprintf(fp, "\n");
    }
}

EventalignSummary summarize_alignment(const SquiggleRead& sr,
                                      uint32_t strand_idx,
                                      const EventAlignmentParameters& params,
                                      const std::vector<EventAlignment>& alignments)
{
    EventalignSummary summary;

    assert(params.alphabet == "");
    const PoreModel* pore_model = params.get_model();
    uint32_t k = pore_model->k;

    size_t prev_ref_pos = std::string::npos;

    // the number of unique reference positions seen in the alignment
    //size_t num_unique_ref_pos = 0;

    for(size_t i = 0; i < alignments.size(); ++i) {

        const EventAlignment& ea = alignments[i];

        summary.num_events += 1;

        // movement information
        size_t ref_move = ea.ref_position - prev_ref_pos;
        if(ref_move == 0) {
            summary.num_stays += 1;
        } else if(i != 0 && ref_move > 1) {
            summary.num_skips += 1;
        } else if(i != 0 && ref_move == 1) {
            summary.num_steps += 1;
        }

        // event information
        summary.sum_duration += sr.get_duration(ea.event_idx, ea.strand_idx);

        if(ea.hmm_state == 'M') {
            uint32_t rank = pore_model->pmalphabet->kmer_rank(ea.model_kmer.c_str(), k);
            double z = z_score(sr, *pore_model, rank, ea.event_idx, ea.strand_idx);
            summary.sum_z_score += z;
        }

        prev_ref_pos = ea.ref_position;
    }

    int nm = bam_aux2i(bam_aux_get(params.record, "NM"));
    summary.alignment_edit_distance = nm;
    if(!alignments.empty()) {
        summary.reference_span = alignments.back().ref_position - alignments.front().ref_position + 1;
    }
    return summary;
}

// Realign the read in event space
void realign_read(EventalignWriter writer,
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

    // load read
    SquiggleRead sr(read_name, read_db, opt::write_samples ? SRF_LOAD_RAW_SAMPLES : 0);

    if(opt::verbose > 1) {
        fprintf(stderr, "Realigning %s [%zu %zu]\n", 
                read_name.c_str(), sr.events[0].size(), sr.events[1].size());
    }
    
    for(int strand_idx = 0; strand_idx < 2; ++strand_idx) {
        
        // Do not align this strand if it was not sequenced
        if(!sr.has_events_for_strand(strand_idx)) {
            continue;
        }

        EventAlignmentParameters params;
        params.sr = &sr;
        params.fai = fai;
        params.hdr = hdr;
        params.record = record;
        params.strand_idx = strand_idx;
        
        params.read_idx = read_idx;
        params.region_start = region_start;
        params.region_end = region_end;

        std::vector<EventAlignment> alignment = align_read_to_ref(params);

        EventalignSummary summary;
        if(writer.summary_fp != NULL) {
            summary = summarize_alignment(sr, strand_idx, params, alignment);
        }

        // write to disk
        #pragma omp critical
        {
            if(opt::output_sam) {
                emit_event_alignment_sam(writer.sam_fp, sr, hdr, record, alignment);
            } else {
                emit_event_alignment_tsv(writer.tsv_fp, sr, strand_idx, params, alignment);
            }

            if(writer.summary_fp != NULL && summary.num_events > 0) {
                assert(params.alphabet == "");
                const PoreModel* pore_model = params.get_model();
                SquiggleScalings& scalings = sr.scalings[strand_idx];
                fprintf(writer.summary_fp, "%zu\t%s\t%s\t", read_idx, read_name.c_str(), sr.fast5_path.c_str());
                fprintf(writer.summary_fp, "%s\t%s\t", pore_model->name.c_str(), strand_idx == 0 ? "template" : "complement");
                fprintf(writer.summary_fp, "%d\t%d\t%d\t%d\t", summary.num_events, summary.num_steps, summary.num_skips, summary.num_stays);
                fprintf(writer.summary_fp, "%.2lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", summary.sum_duration, scalings.shift, scalings.scale, scalings.drift, scalings.var);
            }
        }
    }
}

std::vector<EventAlignment> align_read_to_ref(const EventAlignmentParameters& params)
{
    // Sanity check input parameters
    assert(params.sr != NULL);
    assert(params.fai != NULL);
    assert(params.hdr != NULL);
    assert(params.record != NULL);
    assert(params.strand_idx < NUM_STRANDS);
    assert( (params.region_start == -1 && params.region_end == -1) || (params.region_start <= params.region_end));
    const PoreModel* pore_model = params.get_model();

    std::vector<EventAlignment> alignment_output;

    // Extract the reference subsequence for the entire alignment
    int fetched_len = 0;
    int ref_offset = params.record->core.pos;
    std::string ref_name(params.hdr->target_name[params.record->core.tid]);
    std::string ref_seq = get_reference_region_ts(params.fai, ref_name.c_str(), ref_offset, 
                                                  bam_endpos(params.record), &fetched_len);

    // k from read pore model
    const uint32_t k = params.sr->get_model_k(params.strand_idx);

    // If the reference sequence contains ambiguity codes
    // switch them to the lexicographically lowest base
    ref_seq = pore_model->pmalphabet->disambiguate(ref_seq);
    std::string rc_ref_seq = pore_model->pmalphabet->reverse_complement(ref_seq);

    if(ref_offset == 0)
        return alignment_output;

    // Get the read-to-reference aligned segments
    std::vector<AlignedSegment> aligned_segments = get_aligned_segments(params.record);
    for(size_t segment_idx = 0; segment_idx < aligned_segments.size(); ++segment_idx) {

        AlignedSegment& aligned_pairs = aligned_segments[segment_idx];

        if(params.region_start != -1 && params.region_end != -1) {
            trim_aligned_pairs_to_ref_region(aligned_pairs, params.region_start, params.region_end);
        }

        // Trim the aligned pairs to be within the range of the maximum kmer index
        int max_kmer_idx = params.sr->read_sequence.size() - k;
        trim_aligned_pairs_to_kmer(aligned_pairs, max_kmer_idx);

        if(aligned_pairs.empty())
            return alignment_output;

        bool do_base_rc = bam_is_rev(params.record);
        bool rc_flags[2] = { do_base_rc, !do_base_rc }; // indexed by strand
        const int align_stride = 100; // approximately how many reference bases to align to at once
        const int output_stride = 50; // approximately how many event alignments to output at once

        // get the event range of the read to re-align
        int read_kidx_start = aligned_pairs.front().read_pos;
        int read_kidx_end = aligned_pairs.back().read_pos;
        
        if(do_base_rc) {
            read_kidx_start = params.sr->flip_k_strand(read_kidx_start, k);
            read_kidx_end = params.sr->flip_k_strand(read_kidx_end, k);
        }
        
        assert(read_kidx_start >= 0);
        assert(read_kidx_end >= 0);

        int first_event = params.sr->get_closest_event_to(read_kidx_start, params.strand_idx);
        int last_event = params.sr->get_closest_event_to(read_kidx_end, params.strand_idx);
        bool forward = first_event < last_event;

        int curr_start_event = first_event;
        int curr_start_ref = aligned_pairs.front().ref_pos;
        int curr_pair_idx = 0;

        while( (forward && curr_start_event < last_event) ||
               (!forward && curr_start_event > last_event)) {

            // Get the index of the aligned pair approximately align_stride away
            int end_pair_idx = get_end_pair(aligned_pairs, curr_start_ref + align_stride, curr_pair_idx);
        
            int curr_end_ref = aligned_pairs[end_pair_idx].ref_pos;
            int curr_end_read = aligned_pairs[end_pair_idx].read_pos;

            if(do_base_rc) {
                curr_end_read = params.sr->flip_k_strand(curr_end_read, k);
            }
            assert(curr_end_read >= 0);

            int s = curr_start_ref - ref_offset;
            int l = curr_end_ref - curr_start_ref + 1;

            std::string fwd_subseq = ref_seq.substr(s, l);
            std::string rc_subseq = rc_ref_seq.substr(ref_seq.length() - s - l, l);
            assert(fwd_subseq.length() == rc_subseq.length());

            HMMInputSequence hmm_sequence(fwd_subseq, rc_subseq, pore_model->pmalphabet);
            
            // Require a minimum amount of sequence to align to
            if(hmm_sequence.length() < 2 * k)
                break;

            // Set up HMM input
            HMMInputData input;
            input.read = params.sr;
            input.pore_model = pore_model;
            assert(input.pore_model != NULL);

            input.event_start_idx = curr_start_event;
            input.event_stop_idx = params.sr->get_closest_event_to(curr_end_read, params.strand_idx);
            //printf("[SEGMENT_START] read: %s event start: %zu event end: %zu\n", params.sr->read_name.c_str(), input.event_start_idx, input.event_stop_idx);

            // A limitation of the segment-by-segment alignment is that we can't jump
            // over very large deletions wrt to the reference. The effect of this
            // is that we can get segments that have very few alignable events. We
            // just stop processing them for now
            if(abs((int)input.event_start_idx - (int)input.event_stop_idx) < 2)
                break;

            input.strand = params.strand_idx;
            input.event_stride = input.event_start_idx < input.event_stop_idx ? 1 : -1;
            input.rc = rc_flags[params.strand_idx];

            std::vector<HMMAlignmentState> event_alignment = profile_hmm_align(hmm_sequence, input);
            
            // Output alignment
            size_t num_output = 0;
            size_t event_align_idx = 0;

            // If we aligned to the last event, output everything and stop
            bool last_section = end_pair_idx == (int)aligned_pairs.size() - 1;

            /*
            // Don't allow the segment to end on an E state or else we get alignment
            // artifacts at the segment boundary
            if(!last_section) {
                size_t last_match_index = event_alignment.size() - 1;
                while(event_alignment[last_match_index].state != 'M') {
                    last_match_index -= 1;
                }

                event_alignment.resize(last_match_index + 1);
                if(event_alignment.empty()) {
                    break;
                }
                assert(event_alignment.back().state == 'M');
            }
            */

            int last_event_output = 0;
            int last_ref_kmer_output = 0;

            for(; event_align_idx < event_alignment.size() && 
                  (num_output < output_stride || last_section); event_align_idx++) {

                HMMAlignmentState& as = event_alignment[event_align_idx];
                if(as.state != 'K' && (int)as.event_idx != curr_start_event) {

                    EventAlignment ea;
                    
                    // ref
                    ea.ref_name = ref_name;
                    ea.ref_position = curr_start_ref + as.kmer_idx;
                    ea.ref_kmer = ref_seq.substr(ea.ref_position - ref_offset, k);

                    // event
                    ea.read_idx = params.read_idx;
                    ea.strand_idx = params.strand_idx;
                    ea.event_idx = as.event_idx;
                    ea.rc = input.rc;

                    // hmm
                    ea.hmm_state = as.state;

                    if(ea.hmm_state != 'B') {
                        ea.model_kmer = hmm_sequence.get_kmer(as.kmer_idx, k, input.rc);
                    } else {
                        ea.model_kmer = std::string(k, 'N');
                    }

                    // store
                    alignment_output.push_back(ea);

                    // update
                    last_event_output = as.event_idx;
                    last_ref_kmer_output = curr_start_ref + as.kmer_idx;

                    num_output += 1;
                }
            }

            // Advance the pair iterator to the ref base
            curr_start_event = last_event_output;
            curr_start_ref = last_ref_kmer_output;
            //printf("[SEGMENT_END] read: %s last event output: %zu ref pos: %zu (%s)\n", params.sr->read_name.c_str(), last_event_output, last_ref_kmer_output, ref_seq.substr(last_ref_kmer_output - ref_offset, k).c_str());
            curr_pair_idx = get_end_pair(aligned_pairs, curr_start_ref, curr_pair_idx);

#if EVENTALIGN_TRAIN
            // update training data for read
            params.sr->parameters[params.strand_idx].add_training_from_alignment(hmm_sequence, input, event_alignment);
            global_training[params.strand_idx].add_training_from_alignment(hmm_sequence, input, event_alignment);
#endif

            if(num_output == 0) {
                break;
            }
        } // for realignmentsegment
    } // for bam aligned segment

    return alignment_output;
}

void parse_eventalign_options(int argc, char** argv)
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
            case 'n': opt::print_read_names = true; break;
            case 'f': opt::full_output = true; break;
            case OPT_SAMPLES: opt::write_samples = true; break;
            case 'v': opt::verbose++; break;
            case OPT_MODELS_FOFN: arg >> opt::models_fofn; break;
            case OPT_SCALE_EVENTS: opt::scale_events = true; break;
            case OPT_SUMMARY: arg >> opt::summary_file; break;
            case OPT_SAM: opt::output_sam = true; break;
            case OPT_PROGRESS: opt::progress = true; break;
            case OPT_HELP:
                std::cout << EVENTALIGN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << EVENTALIGN_VERSION_MESSAGE;
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

    if(!opt::models_fofn.empty()) {
        // initialize the model set from the fofn
        PoreModelSet::initialize(opt::models_fofn);
    }

    if (die) 
    {
        std::cout << "\n" << EVENTALIGN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

int eventalign_main(int argc, char** argv)
{
    parse_eventalign_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    ReadDB read_db;
    read_db.load(opt::reads_file);
    
    // Open the BAM and iterate over reads

    // load bam file
    htsFile* bam_fh = sam_open(opt::bam_file.c_str(), "r");
    assert(bam_fh != NULL);

    // load bam index file
    std::string index_filename = opt::bam_file + ".bai";
    hts_idx_t* bam_idx = bam_index_load(index_filename.c_str());
    if(bam_idx == NULL) {
        bam_index_error_exit(opt::bam_file);
    }

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

    // Initialize output
    EventalignWriter writer = { NULL, NULL, NULL };

    if(opt::output_sam) {
        writer.sam_fp = hts_open("-", "w");
        emit_sam_header(writer.sam_fp, hdr);
    } else {
        writer.tsv_fp = stdout;
        emit_tsv_header(writer.tsv_fp);
    }

    if(!opt::summary_file.empty()) {
        writer.summary_fp = fopen(opt::summary_file.c_str(), "w");
        // header
        fprintf(writer.summary_fp, "read_index\tread_name\tfast5_path\tmodel_name\tstrand\tnum_events\t");
        fprintf(writer.summary_fp, "num_steps\tnum_skips\tnum_stays\ttotal_duration\tshift\tscale\tdrift\tvar\n");
    }
    
    // Initialize iteration
    std::vector<bam1_t*> records(opt::batch_size, NULL);
    for(size_t i = 0; i < records.size(); ++i) {
        records[i] = bam_init1();
    }

    int result;
    size_t num_reads_realigned = 0;
    size_t num_records_buffered = 0;
    Progress progress("[eventalign]");

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
                    realign_read(writer, read_db, fai, hdr, record, read_idx, clip_start, clip_end);
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

    if(writer.sam_fp != NULL) {
        hts_close(writer.sam_fp);
    }

    if(writer.summary_fp != NULL) {
        fclose(writer.summary_fp);
    }
    return EXIT_SUCCESS;
}
