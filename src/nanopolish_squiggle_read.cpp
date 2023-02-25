//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_squiggle_read -- Class holding a squiggle (event)
// space nanopore read
//
#include <algorithm>
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_extract.h"
#include "nanopolish_raw_loader.h"
#include "nanopolish_fast5_io.h"

extern "C" {
#include "event_detection.h"
#include "scrappie_common.h"
}

#include <fast5.hpp>

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
SquiggleRead::SquiggleRead(const std::string& name, const ReadDB& read_db, const uint32_t flags) :
    read_name(name),
    nucleotide_type(SRNT_DNA),
    pore_type(PT_UNKNOWN),
    f_p(nullptr)
{

    this->events_per_base[0] = events_per_base[1] = 0.0f;
    this->base_model[0] = this->base_model[1] = NULL;

    if(read_db.get_slow5_mode()){
        slow5_file_t * slow5_file = read_db.get_slow5_file();
        if(!slow5_file){
            fprintf(stderr, "slow5 file is missing");
            exit(EXIT_FAILURE);
        }
        if (!slow5_file->index) {
            fprintf(stderr,"No slow5 index has been loaded\n");
            exit(EXIT_FAILURE);
        }
        slow5_rec_t *rec = NULL;
        int ret = slow5_get(read_name.c_str(), &rec, slow5_file);
        if(ret < 0){
            fprintf(stderr,"Error in when fetching the read\n");
        }

        // metadata
        char* sequencing_kit_c = slow5_hdr_get("sequencing_kit", rec->read_group, slow5_file->header);
        char* experiment_type_c = slow5_hdr_get("experiment_type", rec->read_group, slow5_file->header);

        std::string sequencing_kit = (sequencing_kit_c)?std::string(sequencing_kit_c):"";
        std::string experiment_type = (experiment_type_c)?std::string(experiment_type_c):"";

        char* flowcell_type_c = slow5_hdr_get("flow_cell_product_code", rec->read_group, slow5_file->header);
        // not found, try context_tags + flow_cell_product_code
        if(!flowcell_type_c) {
            flowcell_type_c  = slow5_hdr_get("flowcell_type", rec->read_group, slow5_file->header);
        }
        std::string flowcell_type = (flowcell_type_c)?std::string(flowcell_type_c):"";
        // Flowcell type should always be uppercase
        std::transform(flowcell_type.begin(), flowcell_type.end(), flowcell_type.begin(), ::toupper);

        // Try to detect whether this read is DNA or RNA
        // Fix issue 531: experiment_type in fast5 is "rna" for cDNA kit dcs108
        bool rna_experiment = experiment_type == "rna" || experiment_type == "internal_rna";
        this->nucleotide_type = rna_experiment && sequencing_kit != "sqk-dcs108" ? SRNT_RNA : SRNT_DNA;

        // pre-release R10 uses "cust" flowcell type so we fall back to checking for
        // the pore type in the path. this should be removed eventually once
        // we no longer want to support pre-release R10
        bool cust_or_unknown_fct = flowcell_type.empty() || flowcell_type == "cust-flo-m";
		if( ends_with(sequencing_kit, "114") || flowcell_type == "FLO-MIN110") {
			this->pore_type = PT_R10;
		} else {
			this->pore_type = PT_R9;
		}

        this->read_sequence = read_db.get_read_sequence(read_name);
        load_from_raw_slow5(slow5_file, rec, flags);
        slow5_rec_free(rec);

//        data = Fast5Loader::load_read(slow5_file, name);
    }else{
        this->fast5_path = read_db.get_signal_path(this->read_name);
        g_total_reads += 1;
        if(this->fast5_path == "") {
            g_bad_fast5_file += 1;
            return;
        }

        // Get the read type from the fast5 file
        fast5_file f5_file = fast5_open(fast5_path);
        if(fast5_is_open(f5_file)) {

            std::string sequencing_kit = fast5_get_sequencing_kit(f5_file, this->read_name);
            std::string experiment_type = fast5_get_experiment_type(f5_file, this->read_name);
            std::string flowcell_type = fast5_get_flowcell_type(f5_file, this->read_name);

            // Try to detect whether this read is DNA or RNA
            // Fix issue 531: experiment_type in fast5 is "rna" for cDNA kit dcs108
            bool rna_experiment = experiment_type == "rna" || experiment_type == "internal_rna";
            this->nucleotide_type = rna_experiment && sequencing_kit != "sqk-dcs108" ? SRNT_RNA : SRNT_DNA;

            //fprintf(stderr, "fast5: %s detected flowcell type: %s\n", this->fast5_path.c_str(), flowcell_type.c_str());
            // pre-release R10 uses "cust" flowcell type so we fall back to checking for
            // the pore type in the path. this should be removed eventually once
            // we no longer want to support pre-release R10
            bool cust_or_unknown_fct = flowcell_type.empty() || flowcell_type == "cust-flo-m";
            bool has_r10_path = this->fast5_path.find("r10") != std::string::npos;
            if( ends_with(sequencing_kit, "114") || flowcell_type == "FLO-MIN110") {
                this->pore_type = PT_R10;
            } else {
                this->pore_type = PT_R9;
            }

            this->read_sequence = read_db.get_read_sequence(read_name);
            load_from_raw(f5_file, flags);

            fast5_close(f5_file);

        } else {
            fprintf(stderr, "[warning] fast5 file is unreadable and will be skipped: %s\n", fast5_path.c_str());
            g_bad_fast5_file += 1;
        }
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
static detector_param const event_detection_r10 = {
    .window_length1 = 4,
    .window_length2 = 13,
    .threshold1 = 1.52f,
    .threshold2 = 3.91f,
    .peak_height = 0.17f
};

//
void SquiggleRead::load_from_raw(fast5_file& f5_file, const uint32_t flags)
{
    // File not in db, can't load
    if(this->fast5_path == "" || this->read_sequence == "") {
        return;
    }

    //
    this->read_type = SRT_TEMPLATE;
    std::string strand_str = "template";
    size_t strand_idx = 0;

    // default to R9 parameters
    std::string alphabet = "nucleotide";
    std::string kit = "r9.4_450bps";
    size_t k = 6;
    const detector_param* ed_params = &event_detection_defaults;

    if(this->pore_type == PT_R10) {
        kit = "r10_450bps";
        k = 9;
        //ed_params = &event_detection_r10;
    }

    if(this->nucleotide_type == SRNT_RNA) {
        assert(this->pore_type == PT_R9);
        kit = "r9.4_70bps";
        alphabet = "u_to_t_rna";
        k = 5;
        ed_params = &event_detection_rna;

        std::replace(this->read_sequence.begin(), this->read_sequence.end(), 'U', 'T');
    }


    // Set the base model for this read to either the nucleotide or U->T RNA model
    this->base_model[strand_idx] = PoreModelSet::get_model(kit, alphabet, strand_str, k);
    assert(this->base_model[strand_idx] != NULL);

    // Read the sample rate
    auto channel_params = fast5_get_channel_params(f5_file, this->read_name);
    this->sample_rate = channel_params.sample_rate;

    // Read the actual samples
    raw_table rt = fast5_get_raw_samples(f5_file, this->read_name, channel_params);
    if(rt.n == 0) {
        if(rt.raw != NULL) {
            free(rt.raw);
        }
        return;
    }

    // trim using scrappie's internal method
    // parameters taken directly from scrappie defaults
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
    event_table et = detect_events(rt, *ed_params);
    assert(rt.n > 0);
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
        this->samples.resize(rt.n);
        for(size_t i = 0; i < this->samples.size(); ++i) {
            assert(rt.start + i < rt.n);
            this->samples[i] = rt.raw[rt.start + i];
        }
    }

    // If sequencing RNA, reverse the events to be 5'->3'
    if(this->nucleotide_type == SRNT_RNA) {
        std::reverse(this->events[strand_idx].begin(), this->events[strand_idx].end());
    }

    // clean up scrappie raw and event tables
    assert(rt.raw != NULL);
    assert(et.event != NULL);
    free(rt.raw);
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

void SquiggleRead::load_from_raw_slow5(slow5_file_t* slow5File, slow5_rec_t *rec , const uint32_t flags)
{
    // File not in db, can't load
    if(this->read_sequence == "") {
        return;
    }

    //
    this->read_type = SRT_TEMPLATE;
    std::string strand_str = "template";
    size_t strand_idx = 0;

    // default to R9 parameters
    std::string alphabet = "nucleotide";
    std::string kit = "r9.4_450bps";
    size_t k = 6;
    const detector_param* ed_params = &event_detection_defaults;

    if(this->pore_type == PT_R10) {
        kit = "r10_450bps";
        k = 9;
        //ed_params = &event_detection_r10;
    }

    if(this->nucleotide_type == SRNT_RNA) {
        assert(this->pore_type == PT_R9);
        kit = "r9.4_70bps";
        alphabet = "u_to_t_rna";
        k = 5;
        ed_params = &event_detection_rna;

        std::replace(this->read_sequence.begin(), this->read_sequence.end(), 'U', 'T');
    }


    // Set the base model for this read to either the nucleotide or U->T RNA model
    this->base_model[strand_idx] = PoreModelSet::get_model(kit, alphabet, strand_str, k);
    assert(this->base_model[strand_idx] != NULL);

    // from scrappie

    this->sample_rate = rec->sampling_rate;

    // raw data
    // convert to pA
    float* rawptr = (float*)calloc(rec->len_raw_signal, sizeof(float));
    raw_table rawtbl = { 0, 0, 0, NULL };
    raw_table rt = (raw_table) { rec->len_raw_signal, 0, rec->len_raw_signal, rawptr };
    hsize_t nsample = rec->len_raw_signal;
    float digitisation = rec->digitisation;
    float offset = rec->offset;
    float range = rec->range;
    float raw_unit = range / digitisation;
    for (size_t i = 0; i < nsample; i++) {
        float signal = rec->raw_signal[i];
        rawptr[i] = (signal + offset) * raw_unit;
    }
    if(rt.n == 0) {
        if(rt.raw != NULL) {
            free(rt.raw);
        }
        return;
    }

    // trim using scrappie's internal method
    // parameters taken directly from scrappie defaults
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
    event_table et = detect_events(rt, *ed_params);
    assert(rt.n > 0);
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
        this->samples.resize(rt.n);
        for(size_t i = 0; i < this->samples.size(); ++i) {
            assert(rt.start + i < rt.n);
            this->samples[i] = rt.raw[rt.start + i];
        }
    }

    // If sequencing RNA, reverse the events to be 5'->3'
    if(this->nucleotide_type == SRNT_RNA) {
        std::reverse(this->events[strand_idx].begin(), this->events[strand_idx].end());
    }

    // clean up scrappie raw and event tables
    assert(rt.raw != NULL);
    assert(et.event != NULL);
    free(rt.raw);
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

double median(const std::vector<double>& x)
{
    assert(!x.empty());

    //
    std::vector<double> y = x;
    std::sort(y.begin(), y.end());
    size_t n = y.size();
    return n % 2 == 1 ? y[n / 2] : (y[n / 2 - 1] + y[n / 2]) / 2.0f;
}


RateMetrics SquiggleRead::calculate_rate_metrics(size_t strand_idx) const
{
    RateMetrics rate_metrics;

    // get kmer stats
    size_t k = this->get_base_model(strand_idx)->k;
    size_t num_kmers = this->read_sequence.length() - k + 1;
    size_t num_skips = 0;

    // collect durations, collapsing by k-mer:
    std::vector<double> durations_per_kmer(num_kmers);
    std::vector<double> events_per_kmer(num_kmers);
    for (size_t i = 0; i < this->base_to_event_map.size(); ++i) {
        size_t start_idx = this->base_to_event_map[i].indices[strand_idx].start;
        size_t end_idx = this->base_to_event_map[i].indices[strand_idx].stop;
        // no events for this k-mer
        if (start_idx == -1) {
            num_skips += 1;
            continue;
        }
        assert(start_idx <= end_idx);
        for (size_t j = start_idx; j <= end_idx; ++j) {
            durations_per_kmer[i] += this->get_duration(j, strand_idx);
            events_per_kmer[i] += 1;
        }
    }

    //
    std::sort(durations_per_kmer.begin(), durations_per_kmer.end());
    assert(durations_per_kmer.size() > 0);
    rate_metrics.median_duration = median(durations_per_kmer);

    double stall_threshold = rate_metrics.median_duration * 10;

    double num_any_event = 0;
    double num_extra_event = 0;
    size_t num_stalls = 0;
    double sum_duration;
    double sum_events;
    for(size_t i = 0; i < num_kmers; ++i) {
        sum_duration += durations_per_kmer[i];
        sum_events += events_per_kmer[i];
        num_stalls += durations_per_kmer[i] > stall_threshold;
        num_extra_event += events_per_kmer[i] > 1 ? events_per_kmer[i] - 1 : 0;
        num_any_event += events_per_kmer[i] > 0;
    }
    double mean_duration = sum_duration / num_kmers;
    double mean_events = sum_events / num_kmers;

    // Calculate median duration over segments
    size_t segment_length = 100;
    std::vector<double> segment_mean_duration;
    std::vector<double> segment_mean_events;
    for(size_t i = 0; i < num_kmers; i += segment_length) {
        double segment_sum_duration = 0.0f;
        double segment_sum_events = 0.0f;

        for(size_t j = i; j < num_kmers && j < i + segment_length; ++j) {
            segment_sum_duration += durations_per_kmer[j];
            segment_sum_events += events_per_kmer[j];
        }

        segment_mean_duration.push_back(segment_sum_duration / segment_length);
        segment_mean_events.push_back(segment_sum_events / segment_length);
    }
    double segment_median_duration = median(segment_mean_duration);
    double segment_median_events = median(segment_mean_events);

    // this is our estimator of read rate, currently we use the median duration
    // per k-mer as its more robust to outliers caused by stalls
    double read_rate = 1.0 / rate_metrics.median_duration;
    rate_metrics.skip_frequency = (double)num_skips / num_kmers;
    rate_metrics.stall_frequency = (double)num_stalls / num_kmers;
    //this->extra_event_frequency = (double)num_extra_event / num_any_event;
    rate_metrics.extra_event_frequency = 1 - (1 / (num_extra_event / num_any_event + 1));
    rate_metrics.mean_speed = 1.0 / mean_duration;

    /*
       fprintf(stderr, "[readrate] %s med: [%.5lfs %.1lf] mean: [%.5lfs %.1lf] seg: [%.5lfs %.1f] "
       "ev_mean: %.2f seg med ev: %.2lf, stalls: %zu st f: %.2f skips: %zu skip f: %.4lf ee: %zu ee f: %.4f\n",
       this->read_name.substr(0, 6).c_str(), this->rate_metrics.median_duration, read_rate, mean_duration, 1.0 / mean_duration, segment_median_duration, 1.0 / segment_median_duration,
       mean_events, segment_median_events, num_stalls, (double)num_stalls / num_kmers, num_skips, this->rate_metrics.skip_frequency, num_extra_event, this->rate_metrics.extra_event_frequency);
    */

    /*
    fprintf(stderr, "[read-metadata]\treadid=%s\tmean_read_rate=%.1lf\tmean_events_per_base=%.2lf\tshift=%.2f\tscale=%.2lf\tvar=%.2lf\n",
            this->read_name.c_str(), 1.0 / mean_duration, mean_events, this->scalings[0].shift, this->scalings[0].scale, this->scalings[0].var);
    */
    return rate_metrics;
}

       
SNRMetrics SquiggleRead::calculate_snr_metrics(size_t strand_idx) const
{
    SNRMetrics out = { 0.0f, 0.0f };
    if(this->events[strand_idx].size() < 100) {
        return out;
    }

    // SNR estimates (from CW ONT)
    std::vector<double> event_means;
    std::vector<double> event_sds;

    for(size_t i = 0; i < this->events[strand_idx].size(); ++i) {
        event_means.push_back(this->events[strand_idx][i].mean);
        event_sds.push_back(this->events[strand_idx][i].stdv);
    }

    std::sort(event_means.begin(), event_means.end());
    std::sort(event_sds.begin(), event_sds.end());

    int idx10p = event_means.size() * 0.1;
    int idx90p = event_means.size() * 0.9;
    out.current_range = event_means[idx90p] - event_means[idx10p];
    
    int idx50p = event_sds.size() * 0.5;
    out.median_sd = event_sds[idx50p];
    return out;
}
