//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_squiggle_read -- Class holding a squiggle (event)
// space nanopore read
//
#include <algorithm>
#include <limits>
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_raw_loader.h"
#include "nanopolish_fast5_io.h"
#include "nanopolish_fast5_loader.h"

extern "C" {
#include "event_detection.h"
#include "scrappie_common.h"
}

//#define DEBUG_MODEL_SELECTION 1
//#define DEBUG_RECONSTRUCTION 1

// Track the number of skipped reads to warn the use at the end of the run
// Workaround for albacore issues.  Temporary, I hope
int g_total_reads = 0;
int g_unparseable_reads = 0;
int g_qc_fail_reads = 0;
int g_failed_calibration_reads = 0;
int g_failed_alignment_reads = 0;
int g_bad_fast5_file = 0;

const double MIN_CALIBRATION_VAR = 2.5;

void SquiggleScalings::set4(double _shift,
                            double _scale,
                            double _drift,
                            double _var)
{
    set6(_shift, _scale, _drift, _var, 1.0, 1.0);
}

void SquiggleScalings::set6(double _shift,
                            double _scale,
                            double _drift,
                            double _var,
                            double _scale_sd,
                            double _var_sd)
{
    // direct
    shift = _shift;
    scale = _scale;
    drift = _drift;
    var = _var;
    scale_sd = _scale_sd;
    var_sd = _var_sd;

    // derived
    log_var = log(var);
    scaled_var = var / scale;
    log_scaled_var = log(scaled_var);
}

//
SquiggleRead::SquiggleRead(const std::string& name, const ReadDB& read_db, const uint32_t flags)
{
    Fast5Data data;
    std::string sequence = read_db.get_read_sequence(name);
    if(sequence.empty()){
        fprintf(stderr,"[warning] sequence of read %s is empty\n", name.c_str());
        return;
    }

    if(read_db.get_slow5_mode()) {
        slow5_file_t* slow5_file = read_db.get_slow5_file();
        if(!slow5_file) {
            fprintf(stderr, "slow5 file is missing");
            exit(EXIT_FAILURE);
        }
        if(!slow5_file->index) {
            fprintf(stderr,"No slow5 index has been loaded\n");
            exit(EXIT_FAILURE);
        }
        data = Fast5Loader::load_read(slow5_file, name);
    } else {
        this->fast5_path = read_db.get_signal_path(name);
        if(this->fast5_path == "") {
            g_bad_fast5_file += 1;
            return;
        }
        data = Fast5Loader::load_read(fast5_path, name);
        if(!data.is_valid) {
            fprintf(stderr, "[warning] fast5 file is unreadable and will be skipped: %s\n", fast5_path.c_str());
            g_bad_fast5_file += 1;
        }
    }
    if(data.is_valid) {
        init(sequence, data, flags);
    }
    if(!this->events[0].empty()) {
        assert(this->base_model[0] != NULL);
    }
    free(data.rt.raw);
    data.rt.raw = NULL;
}

SquiggleRead::SquiggleRead(const ReadDB& read_db, const Fast5Data& data, const uint32_t flags)
{
    init(read_db.get_read_sequence(data.read_name), data, flags);
}

SquiggleRead::SquiggleRead(const std::string& sequence, const Fast5Data& data, const uint32_t flags)
{
    init(sequence, data, flags);
}

//
void SquiggleRead::init(const std::string& read_sequence, const Fast5Data& data, const uint32_t flags)
{
    this->nucleotide_type = SRNT_DNA;
    this->pore_type = PORETYPE_UNKNOWN;

    this->events_per_base[0] = events_per_base[1] = 0.0f;
    this->base_model[0] = this->base_model[1] = NULL;
    g_total_reads += 1;

    this->read_name = data.read_name;
    this->read_sequence = read_sequence;

    // sometimes the basecaller will emit very short sequences, which causes problems
    // also there can be rare issues with the signal in the fast5 and we want to skip
    // such reads
    if(this->read_sequence.length() > 20 && data.is_valid && data.rt.n > 0) {
        load_from_raw(data, flags);
    } else {
        g_bad_fast5_file += 1;
    }

    if(!this->events[0].empty()) {
        assert(this->base_model[0] != NULL);
    }
}

SquiggleRead::~SquiggleRead()
{

}

// helper for get_closest_event_to
int SquiggleRead::get_next_event(int start, int stop, int stride, uint32_t strand) const
{
    while(start != stop) {

        int ei = base_to_event_map[start].indices[strand].start;
        if(ei != -1)
            return ei;
        start += stride;
    }
    return -1;
}

//
int SquiggleRead::get_closest_event_to(int k_idx, uint32_t strand) const
{
    int stop_before = std::max(0, k_idx - 1000);
    int stop_after = std::min(k_idx + 1000, (int)base_to_event_map.size() - 1);

    int event_before = get_next_event(k_idx, stop_before, -1, strand);
    int event_after = get_next_event(k_idx, stop_after, 1, strand);

    // TODO: better selection of "best" event to return
    if(event_before == -1)
        return event_after;
    return event_before;
}

//
void SquiggleRead::load_from_raw(const Fast5Data& fast5_data, const uint32_t flags)
{

    // Try to detect whether this read is DNA or RNA
    // Fix issue 531: experiment_type in fast5 is "rna" for cDNA kit dcs108
    bool rna_experiment = fast5_data.experiment_type == "rna" || fast5_data.experiment_type == "internal_rna";
    this->nucleotide_type = rna_experiment && fast5_data.sequencing_kit != "sqk-dcs108" ? SRNT_RNA : SRNT_DNA;

    // Hardcoded parameters, for now we can only do template with the main R9.4 model
    size_t strand_idx = 0;
    std::string alphabet = "nucleotide";
    std::string kit = "r9.4_450bps";
    std::string strand_str = "template";
    size_t k = 6;

    const detector_param* ed_params = &event_detection_defaults;

    if(this->nucleotide_type == SRNT_RNA) {
        kit = "r9.4_70bps";
        alphabet = "u_to_t_rna";
        k = 5;
        ed_params = &event_detection_rna;

        std::replace(this->read_sequence.begin(), this->read_sequence.end(), 'U', 'T');
    }

    this->read_type = SRT_TEMPLATE;
    this->pore_type = PORETYPE_R9;

    // Set the base model for this read to either the nucleotide or U->T RNA model
    this->base_model[strand_idx] = PoreModelSet::get_model(kit, alphabet, strand_str, k);
    assert(this->base_model[strand_idx] != NULL);

    // Read the sample rate
    this->sample_rate = fast5_data.channel_params.sample_rate;
    this->channel_id = fast5_data.channel_params.channel_id;
    this->sample_start_time = fast5_data.start_time;

    // trim raw using scrappie's internal method
    // parameters taken directly from scrappie defaults
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(fast5_data.rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
    event_table et = detect_events(fast5_data.rt, *ed_params);
    assert(et.n > 0);

    //
    this->scalings[strand_idx] = estimate_scalings_using_mom(this->read_sequence,
                                                             *this->base_model[strand_idx],
                                                             et);

    // copy events into nanopolish's format
    this->events[strand_idx].resize(et.n);
    double start_time = 0;
    for(size_t i = 0; i < et.n; ++i) {
        float length_in_seconds = et.event[i].length / this->sample_rate;
        this->events[strand_idx][i] = { et.event[i].mean, et.event[i].stdv, start_time, length_in_seconds, logf(et.event[i].stdv) };
        start_time += length_in_seconds;
    }

    if(flags & SRF_LOAD_RAW_SAMPLES) {
        this->sample_start_time = 0;
        this->samples.resize(fast5_data.rt.n);
        for(size_t i = 0; i < this->samples.size(); ++i) {
            assert(fast5_data.rt.start + i < fast5_data.rt.n);
            this->samples[i] = fast5_data.rt.raw[fast5_data.rt.start + i];
        }
    }

    // If sequencing RNA, reverse the events to be 5'->3'
    if(this->nucleotide_type == SRNT_RNA) {
        std::reverse(this->events[strand_idx].begin(), this->events[strand_idx].end());
    }

    // clean up event tables
    assert(et.event != NULL);
    free(et.event);

    // align events to the basecalled read
    std::vector<AlignedPair> event_alignment = adaptive_banded_simple_event_align(*this, *this->base_model[strand_idx], read_sequence);

    // transform alignment into the base-to-event map
    if(event_alignment.size() > 0) {

        // create base-to-event map
        size_t n_kmers = read_sequence.size() - this->get_model_k(strand_idx) + 1;
        this->base_to_event_map.clear();
        this->base_to_event_map.resize(n_kmers);

        size_t max_event = 0;
        size_t min_event = std::numeric_limits<size_t>::max();

        size_t prev_event_idx = -1;
        for(size_t i = 0; i < event_alignment.size(); ++i) {

            size_t k_idx = event_alignment[i].ref_pos;
            size_t event_idx = event_alignment[i].read_pos;
            IndexPair& elem = this->base_to_event_map[k_idx].indices[strand_idx];
            if(event_idx != prev_event_idx) {
                if(elem.start == -1) {
                    elem.start = event_idx;
                }
                elem.stop = event_idx;
            }

            max_event = std::max(max_event, event_idx);
            min_event = std::min(min_event, event_idx);
            prev_event_idx = event_idx;
        }

        events_per_base[strand_idx] = (double)(max_event - min_event) / n_kmers;

        // prepare data structures for the final calibration
        std::vector<EventAlignment> alignment =
            get_eventalignment_for_1d_basecalls(read_sequence, alphabet, this->base_to_event_map, this->base_model[strand_idx]->k, strand_idx, 0);

        // run recalibration to get the best set of scaling parameters and the residual
        // between the (scaled) event levels and the model.
        // internally this function will set shift/scale/etc of the pore model
        bool calibrated = recalibrate_model(*this, *this->base_model[strand_idx], strand_idx, alignment, true, false);

#ifdef DEBUG_MODEL_SELECTION
        fprintf(stderr, "[calibration] read: %s events: %zu"
                         " scale: %.2lf shift: %.2lf drift: %.5lf var: %.2lf\n",
                                read_name.substr(0, 6).c_str(), this->events[strand_idx].size(), this->scalings[strand_idx].scale,
                                this->scalings[strand_idx].shift, this->scalings[strand_idx].drift, this->scalings[strand_idx].var);
#endif

        // QC calibration
        if(!calibrated || this->scalings[strand_idx].var > MIN_CALIBRATION_VAR) {
            events[strand_idx].clear();
            g_failed_calibration_reads += 1;
        }
    } else {
        // Could not align, fail this read
        this->events[strand_idx].clear();
        this->events_per_base[strand_idx] = 0.0f;
        g_failed_alignment_reads += 1;
    }

    // Filter poor quality reads that have too many "stays"
    if(!this->events[strand_idx].empty() && this->events_per_base[strand_idx] > 5.0) {
        g_qc_fail_reads += 1;
        events[0].clear();
        events[1].clear();
    }
}

// Return a vector of eventalignments for the events that made up the basecalls in the read
std::vector<EventAlignment> SquiggleRead::get_eventalignment_for_1d_basecalls(const std::string& read_sequence_1d,
                                                                              const std::string& alphabet_name,
                                                                              const std::vector<EventRangeForBase>& base_to_event_map_1d,
                                                                              const size_t k,
                                                                              const size_t strand_idx,
                                                                              const int shift_offset) const
{
    std::vector<EventAlignment> alignment;

    const Alphabet* alphabet = get_alphabet_by_name(alphabet_name);
    size_t n_kmers = read_sequence_1d.size() - k + 1;
    size_t prev_kmer_rank = -1;

    for(int ki = 0; ki < n_kmers; ++ki) {
        IndexPair event_range_for_kmer = base_to_event_map_1d[ki].indices[strand_idx];

        // skip kmers without events
        if(event_range_for_kmer.start == -1)
            continue;

        // skip k-mers that cannot be shifted to a valid position
        if(ki + shift_offset < 0 || ki + shift_offset >= n_kmers) {
            continue;
        }

        for(size_t event_idx = event_range_for_kmer.start;
            event_idx <= event_range_for_kmer.stop; event_idx++)
        {
            assert(event_idx < this->events[strand_idx].size());

            // since we use the 1D read seqence here we never have to reverse complement
            std::string kmer = read_sequence_1d.substr(ki + shift_offset, k);
            size_t kmer_rank = alphabet->kmer_rank(kmer.c_str(), k);

            EventAlignment ea;
            // ref data
            //ea.ref_name = "read";
            ea.read_idx = -1; // not needed
            ea.ref_kmer = kmer;
            ea.ref_position = ki;
            ea.strand_idx = strand_idx;
            ea.event_idx = event_idx;
            ea.rc = false;
            ea.model_kmer = kmer;
            ea.hmm_state = prev_kmer_rank != kmer_rank ? 'M' : 'E';
            alignment.push_back(ea);
            prev_kmer_rank = kmer_rank;
        }
    }

    return alignment;
}

size_t SquiggleRead::get_sample_index_at_time(size_t sample_time) const
{
    return sample_time - this->sample_start_time;
}

//
std::vector<float> SquiggleRead::get_scaled_samples_for_event(size_t strand_idx, size_t event_idx) const
{
    std::pair<size_t, size_t> sample_range = get_event_sample_idx(strand_idx, event_idx);

    std::vector<float> out;
    for(size_t i = sample_range.first; i < sample_range.second; ++i) {
        double curr_sample_time = (this->sample_start_time + i) / this->sample_rate;
        //fprintf(stderr, "event_start: %.5lf sample start: %.5lf curr: %.5lf rate: %.2lf\n", event_start_time, this->sample_start_time / this->sample_rate, curr_sample_time, this->sample_rate);
        double s = this->samples[i];
        // apply scaling corrections
        double scaled_s = s - this->scalings[strand_idx].shift;
        assert(curr_sample_time >= (this->sample_start_time / this->sample_rate));
        scaled_s -= (curr_sample_time - (this->sample_start_time / this->sample_rate)) * this->scalings[strand_idx].drift;
        scaled_s /= this->scalings[strand_idx].scale;
        out.push_back(scaled_s);
    }
    return out;
}

// return a pair of value corresponding to the start and end index of a given index on the signal
std::pair<size_t, size_t> SquiggleRead::get_event_sample_idx(size_t strand_idx, size_t event_idx) const
{
    double event_start_time = this->events[strand_idx][event_idx].start_time;
    double event_duration = this->events[strand_idx][event_idx].duration;

    size_t start_idx = this->get_sample_index_at_time(event_start_time * this->sample_rate);
    size_t end_idx = this->get_sample_index_at_time((event_start_time + event_duration) * this->sample_rate);

    return std::make_pair(start_idx, end_idx);
}
