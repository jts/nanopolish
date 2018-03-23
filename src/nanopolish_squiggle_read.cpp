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
    this->fast5_path = read_db.get_signal_path(this->read_name);
    g_total_reads += 1;
    if(this->fast5_path == "") {
        g_bad_fast5_file += 1;
        return;
    }


    // Get the read type from the fast5 file
    hid_t hdf5_file = fast5_open(fast5_path);
    if(hdf5_file >= 0) {

        //fprintf(stderr, "file: %s\n", fast5_path.c_str());
        std::string experiment_type = fast5_get_experiment_type(hdf5_file);
        //fprintf(stderr, "type: %s\n", experiment_type.c_str());

        // Try to detect whether this read is DNA or RNA
        this->nucleotide_type = experiment_type == "rna" ? SRNT_RNA : SRNT_DNA;

        // Did this read come from nanopolish extract?
        bool is_event_read = is_extract_read_name(this->read_name);

        // Use the legacy loader for DNA reads from extract, otherwise use the new loader
        if(this->nucleotide_type == SRNT_DNA && is_event_read) {
            try {
                #pragma omp critical(sr_load_fast5)
                {
                    this->f_p = new fast5::File(fast5_path);
                    assert(this->f_p->is_open());
                    load_from_events(flags);
                }
            } catch(hdf5_tools::Exception e) {
                fprintf(stderr, "[warning] fast5 file is unreadable and will be skipped: %s\n", fast5_path.c_str());
                g_bad_fast5_file += 1;
            }

            delete this->f_p;
            this->f_p = nullptr;
        } else {
            this->read_sequence = read_db.get_read_sequence(read_name);
            load_from_raw(hdf5_file, flags);
        }

        fast5_close(hdf5_file);

    } else {
        fprintf(stderr, "[warning] fast5 file is unreadable and will be skipped: %s\n", fast5_path.c_str());
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
void SquiggleRead::load_from_events(const uint32_t flags)
{
    assert(this->nucleotide_type != SRNT_RNA);

    assert(f_p->is_open());
    detect_pore_type();
    detect_basecall_group();
    assert(not basecall_group.empty());

    read_sequence = f_p->get_basecall_seq(read_type, basecall_group);

    // Load PoreModel for both strands
    std::vector<EventRangeForBase> event_maps_1d[NUM_STRANDS];
    std::string read_sequences_1d[NUM_STRANDS];

    for (size_t si = 0; si < 2; ++si) {

        // Do we want to load this strand?
        if(! (read_type == SRT_2D || read_type == si) ) {
            continue;
        }

        // Load the events for this strand
        auto f5_events = f_p->get_basecall_events(si, basecall_group);

        // copy events
        events[si].resize(f5_events.size());
        std::vector<double> p_model_states;

        for(size_t ei = 0; ei < f5_events.size(); ++ei) {
            auto const & f5_event = f5_events[ei];

            events[si][ei] = { static_cast<float>(f5_event.mean),
                               static_cast<float>(f5_event.stdv),
                               f5_event.start,
                               static_cast<float>(f5_event.length),
                               static_cast<float>(log(f5_event.stdv))
                             };
            assert(f5_event.p_model_state >= 0.0 && f5_event.p_model_state <= 1.0);
            p_model_states.push_back(f5_event.p_model_state);
        }

        // we need the 1D event map and sequence to calculate calibration parameters
        // these will be copied into the member fields later if this is a 1D read,
        // or discarded if this is a 2D read

        // NB we use event_group in this call rather than basecall_group as we want the 1D basecalls that match the events
        read_sequences_1d[si] = f_p->get_basecall_seq(si == 0 ? SRT_TEMPLATE : SRT_COMPLEMENT,
                                                      f_p->get_basecall_1d_group(basecall_group));
        event_maps_1d[si] = build_event_map_1d(read_sequences_1d[si], si, f5_events, 5);

        // Constructing the event map can fail due to an albacore bug
        // in this case, we have to set this strand to be invalid
        if(!event_maps_1d[si].empty()) {
            // run version-specific load
            if(pore_type == PT_R7) {
                _load_R7(si);
            } else {
                _load_R9(si, read_sequences_1d[si], event_maps_1d[si], p_model_states, flags);
            }
        } else {
            events[si].clear();
        }
    }

    // Build the map from k-mers of the read sequence to events
    if(read_type == SRT_2D) {
        if(pore_type == PT_R9) {
            build_event_map_2d_r9();
        } else {
            assert(pore_type == PT_R7);
            build_event_map_2d_r7();
        }
    } else {
        assert(read_type < NUM_STRANDS);
        this->base_to_event_map.swap(event_maps_1d[read_type]);
    }

    // Load raw samples if requested
    if(flags & SRF_LOAD_RAW_SAMPLES) {

        auto& sample_read_names = f_p->get_raw_samples_read_name_list();
        if(sample_read_names.empty()) {
            fprintf(stderr, "Error, no raw samples found\n");
            exit(EXIT_FAILURE);
        }

        // we assume the first raw sample read is the one we're after
        std::string sample_read_name = sample_read_names.front();

        samples = f_p->get_raw_samples(sample_read_name);
        sample_start_time = f_p->get_raw_samples_params(sample_read_name).start_time;

        // retreive parameters
        auto channel_params = f_p->get_channel_id_params();
        sample_rate = channel_params.sampling_rate;
    }

    // Filter poor quality reads that have too many "stays"
    if(!events[0].empty() && events_per_base[0] > 5.0) {
        g_qc_fail_reads += 1;
        events[0].clear();
        events[1].clear();
    }
}

//
void SquiggleRead::load_from_raw(hid_t hdf5_file, const uint32_t flags)
{
    // File not in db, can't load
    if(this->fast5_path == "" || this->read_sequence == "") {
        return;
    }

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
    this->pore_type = PT_R9;

    // Set the base model for this read to either the nucleotide or U->T RNA model
    this->base_model[strand_idx] = PoreModelSet::get_model(kit, alphabet, strand_str, k);
    assert(this->base_model[strand_idx] != NULL);

    // Read the sample rate
    auto channel_params = fast5_get_channel_params(hdf5_file);
    this->sample_rate = channel_params.sample_rate;

    // Read the actual samples
    raw_table rt = fast5_get_raw_samples(hdf5_file, channel_params);

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

    // If sequencing RNA, reverse the events to be 3'->5'
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

void SquiggleRead::_load_R7(uint32_t si)
{
    assert(f_p and f_p->is_open());
    assert(not basecall_group.empty());

    // Load the pore model for this strand
    PoreModel tmp_pm = PoreModel( f_p, si, basecall_group );
    if(!PoreModelSet::has_model(tmp_pm)) {
        this->base_model[si] = PoreModelSet::add_model(tmp_pm);
    }

    // Read scaling parameters
    auto params = f_p->get_basecall_model_params(si, basecall_group);
    this->scalings[si].set6(params.shift, params.scale, params.drift, params.var, params.scale_sd, params.var_sd);

    // initialize transition parameters
    parameters[si].initialize(this->base_model[si]->metadata);
}

void SquiggleRead::_load_R9(uint32_t si,
                            const std::string& read_sequence_1d,
                            const std::vector<EventRangeForBase>& event_map_1d,
                            const std::vector<double>& p_model_states,
                            const uint32_t flags)
{
    size_t calibration_k = 5;
    size_t final_model_k = 6;

    assert(f_p and f_p->is_open());
    
    // The k-mer label semantics differ between basecaller versions and models
    // We use a "label_shift" parameter to determine how to line up the labels
    // with events so that we can recalibrate the models.
    fast5::Attr_Map basecall_attributes = f_p->get_basecall_params(basecall_group);
    std::string basecaller_name = basecall_attributes["name"];

    bool is_albacore = false;
    bool is_albacore_1_or_later = false;
    if(basecaller_name.find("Albacore") != -1) {
        is_albacore = true;
        SemVer ver = parse_semver_string(basecall_attributes["version"]);
        if(ver.major >= 1) {
            is_albacore_1_or_later = true;
        }
    }

    // Parse kit name and label_shift from the model type encoded in the fast5
    std::string kit = "";
    int label_shift = 0;
    if( (flags & SRF_NO_MODEL) == 0) {

        auto config = f_p->get_basecall_config(basecall_group);
        std::string mt = "";

        if(is_albacore) {

            // Deal with albacore 1.1.0 changing the model convention
            mt = config["basecall_1d/model"];
            if(mt == "") {
                mt = config["basecall_1d/template_model"];
            }

            // remove prefix/suffix
            auto fields = split(mt, '_');
            if(! (fields.size() == 4 || (fields.size() == 5 && fields[4] == "nontransducer.jsn")) ) {
                fprintf(stderr, "error: could not parse model string: %s\n", mt.c_str());
                exit(EXIT_FAILURE);
            }
            mt = fields[1] + "_" + fields[2];
        } else {
            mt = config["general/model_type"];
        }

        if(mt == "") {
            fprintf(stderr, "Error: The basecalling model could not be detected from the fast5 file: %s\n", this->fast5_path.c_str());
            fprintf(stderr, "Error: Please re-run basecalling using albacore.\n");
            exit(1);
        }

        kit = "r9.4_450bps";
        // all 250bps data should use this model (according to ONT see
        // https://github.com/nanoporetech/kmer_models/issues/3)
        if(mt == "r9_250bps_nn" || mt == "r9_250bps" || mt == "r94_250bps" || mt == "r94_250bps_nn" || mt == "r9.4_250bps") {
            label_shift = 0;
            kit = "r9_250bps";
        } else if(mt == "r94_450bps" || mt == "r9_450bps" || mt == "r9.4_450bps" || mt == "r9.5_450bps") {
            label_shift = is_albacore_1_or_later ? -1 : 0;
            kit = "r9.4_450bps";
        } else {
            fprintf(stderr, "Unknown model type string: %s, please report on github.\n", mt.c_str());
            exit(1);
        }
    }

    std::vector<EventAlignment> alignment =
        get_eventalignment_for_1d_basecalls(read_sequence_1d, "nucleotide", event_map_1d, calibration_k, si, label_shift);

    // JTS Hack: blacklist bad k-mer and filter out events with low p_model_state
    double keep_fraction = 0.75;
    std::vector<double> sorted_p_model_states = p_model_states;
    std::sort(sorted_p_model_states.begin(), sorted_p_model_states.end());
    double p_model_state_threshold = sorted_p_model_states[sorted_p_model_states.size() * (1 - keep_fraction)];

    std::string blacklist_kmer = "CCTAG";
    std::vector<EventAlignment> filtered;
    filtered.reserve(alignment.size());

    assert(p_model_states.size() == events[si].size());

    // This vector tracks the number of events observed (by the basecaller) for each kmer
    std::vector<size_t> event_counts;
    event_counts.reserve(read_sequence_1d.length());
    std::string prev_kmer = "";

    for(const auto& ea : alignment) {
        if((!ea.rc && ea.ref_kmer == blacklist_kmer) ||
           (ea.rc && ea.ref_kmer == gDNAAlphabet.reverse_complement(blacklist_kmer)))
        {
            continue;
        }

        if(ea.ref_kmer != prev_kmer) {
            prev_kmer = ea.ref_kmer;
            event_counts.push_back(1);
        } else {
            assert(!event_counts.empty());
            event_counts.back() += 1;
        }

        if(p_model_states[ea.event_idx] < p_model_state_threshold)
            continue;

        filtered.push_back(ea);
    }


    // Estimate the events-per-base sequencing rate
    // This is used for QC and to parameterize the HMM
    // We trim off the lowest and highest trim_frac event counts
    // to avoid the extremely long stays that occasionally occur
    std::sort(event_counts.begin(), event_counts.end());
    double trim_frac = 0.05;
    size_t trim_start_idx = trim_frac * (float)event_counts.size();
    size_t trim_end_idx = (1 - trim_frac) * (float)event_counts.size();
    size_t event_count_sum = 0;
    for(size_t i = trim_start_idx; i < trim_end_idx; ++i) {
        event_count_sum += event_counts[i];
    }

    // Estimate sequencing rate
    double total_duration = get_time(events[0].size() - 1, 0);
    double rate = read_sequence_1d.size() / total_duration;

    events_per_base[si] = (float)event_count_sum / (trim_end_idx - trim_start_idx);
    //fprintf(stderr, "events per base: %.2lf rate: %.2lf\n", events_per_base[si], rate);

    // Load the pore model (if requested) and calibrate it
    if( (flags & SRF_NO_MODEL) == 0) {

        std::string alphabet = "nucleotide"; // always calibrate with the nucleotide alphabet

        // For the template strad we only have one candidate model
        // For complement we need to select between the two possible models
        std::vector<const PoreModel*> candidate_models;
        if(si == 0) {
            candidate_models.push_back(PoreModelSet::get_model(kit, alphabet, "template", calibration_k));
        } else {
            for(const std::string& cmn : { "complement.pop1", "complement.pop2" } ) {
                if(PoreModelSet::has_model(kit, alphabet, cmn, calibration_k)) {
                   candidate_models.push_back(PoreModelSet::get_model(kit, alphabet, cmn, calibration_k));
                }
            }
        }

        const PoreModel* best_model;
        SquiggleScalings best_scalings;
        double best_model_var = INFINITY;

        for(size_t model_idx = 0; model_idx < candidate_models.size(); model_idx++) {

            const PoreModel* curr_model = candidate_models[model_idx];

            // run recalibration to get the best set of scaling parameters and the residual
            // between the (scaled) event levels and the model
            bool calibrated = recalibrate_model(*this, *curr_model, si, filtered, true, false);
            if(calibrated) {
                if(this->scalings[si].var < best_model_var) {
                    best_model_var = this->scalings[si].var;
                    best_model = curr_model;
                    best_scalings = this->scalings[si];
                }
            }

#ifdef DEBUG_MODEL_SELECTION
            fprintf(stderr, "[calibration] read: %s strand: %zu model_idx: %zu "
                             "scale: %.2lf shift: %.2lf drift: %.5lf var: %.2lf\n",
                                    read_name.substr(0, 6).c_str(), si, model_idx, scalings[si].scale,
                                    scalings[si].shift, scalings[si].drift, scalings[si].var);
#endif
        }

        if(best_model_var < MIN_CALIBRATION_VAR) {
#ifdef DEBUG_MODEL_SELECTION
            fprintf(stderr, "[calibration] selected model with var %.4lf\n", best_model_var);
#endif
            // Replace model
            const PoreModel* final_model = PoreModelSet::get_model(kit,
                                                                   alphabet,
                                                                   best_model->metadata.get_strand_model_name(),
                                                                   final_model_k);
 
            this->base_model[si] = final_model;
            this->scalings[si] = best_scalings;

            // initialize transition parameters
            parameters[si].initialize(this->base_model[si]->metadata);
        } else {
            g_failed_calibration_reads += 1;
            // could not find a model for this strand, discard it
            events[si].clear();
        }
    }
}

inline size_t search_for_event_kmer(const std::string sequence,
                                    size_t start,
                                    size_t num_seq_kmers,
                                    size_t k,
                                    const std::string& event_kmer,
                                    size_t max_steps)
{
    for(size_t i = 0; i + start < num_seq_kmers && i <= max_steps; ++i) {
        if(sequence.compare(start + i, k, event_kmer, 0, k) == 0) {
            return i;
        }
    }
    return max_steps + 1;
}


size_t assign_events_for_segment(const std::string& read_sequence_1d,
                               const std::string& segment,
                               const uint32_t k,
                               size_t previous_segment_kmer_end,
                               size_t start_event_idx,
                               const std::vector<fast5::Basecall_Event>& f5_events,
                               std::vector<EventRangeForBase>& out_event_map)
{
    // Find the start position of this segment in the basecalled read
    // To avoid searching the entire read sequence we start a few bp
    // upstream of the end of the last segment
    const size_t backtrack_distance = 20;
    size_t start_hint = 0;
    if(previous_segment_kmer_end >= backtrack_distance) {
        start_hint = previous_segment_kmer_end - backtrack_distance;
    }
    size_t start_kmer_idx = read_sequence_1d.find(segment, start_hint);
    size_t end_kmer_idx = start_kmer_idx + segment.length() - k + 1;
    fprintf(stderr, "[segment sequence length=%zu start=%zu end=%zu] %s\n", segment.length(), start_kmer_idx, end_kmer_idx, segment.c_str());
    assert(start_kmer_idx < read_sequence_1d.size());
    return end_kmer_idx;
}


void segment_by_moves(const std::string& read_sequence_1d,
                      const std::vector<fast5::Basecall_Event>& f5_events,
                      std::vector<EventRangeForBase>& out_event_map)
{
    assert(!f5_events.empty());
    std::string segment = "";
    size_t sum_move = 0;
    size_t segment_event_start = 0;
    size_t previous_segment_kmer_end = 0;
    const uint32_t k = array2str(f5_events.front().model_state).length();
    size_t prefix_skip_length = 0;
    std::string previous_kmer = "";

    for(size_t curr_event_idx = 0; curr_event_idx < f5_events.size(); ++curr_event_idx) {

        std::string event_kmer = array2str(f5_events[curr_event_idx].model_state);
        size_t k = event_kmer.length();
        if(prefix_skip_length > 0) {
            // Skip over this k-mer, decrease counter if k-mer is different than previous
            if(event_kmer != previous_kmer) {
                previous_kmer = event_kmer;
                prefix_skip_length -= 1;
            }
            continue;
        }

        // Determine overlap length with previous k-mer
        // If this is the first k-mer in the segment we set the overlap length to be zero
        // and ignore the move field
        uint32_t move;
        if(segment.empty()) {
            segment_event_start = curr_event_idx;
            move = k;
        } else {
            move = f5_events[curr_event_idx].move;
        }

        int overlap = k - move;

        // Compute the overlap and overhang string
        std::string overlap_str = event_kmer.substr(0, overlap);
        std::string overhang_str = event_kmer.substr(overlap);
        sum_move += move;
        uint32_t display_suffix_length = std::min(k, segment.length());

        size_t read_kmer_idx = previous_segment_kmer_end + sum_move - k;
        std::string read_kmer = read_sequence_1d.substr(read_kmer_idx, k);
        //fprintf(stderr, "[%zu] seq suffix: %s move: %zu overlap: %s overhang: %s read: %s\n", read_kmer_idx, segment.substr(segment.length() - display_suffix_length).c_str(), move, overlap_str.c_str(), overhang_str.c_str(), read_kmer.c_str());
        // Check if the overlap between k-mers is valid. If it is not then we end the segment
        if(overlap_str == segment.substr(segment.length() - overlap)) {
            if(read_kmer != event_kmer) {
                //fprintf(stderr, "read kmer does not match event k-mer despite good overlap\n");
            }
            segment.append(overhang_str);
        } else {

            // up to k k-mers may be incorrect in the segment, trim them off
            std::string trimmed_segment = segment.substr(0, segment.length() - k);
            previous_segment_kmer_end = assign_events_for_segment(read_sequence_1d, trimmed_segment, k, previous_segment_kmer_end, segment_event_start, f5_events, out_event_map);

            // reset variables
            segment = "";
            segment_event_start = curr_event_idx;
            prefix_skip_length = k;
            sum_move = 0;
        }
    }

    // Handle the last segment
    assign_events_for_segment(read_sequence_1d, segment, k, previous_segment_kmer_end, segment_event_start, f5_events, out_event_map);
}

std::vector<EventRangeForBase> SquiggleRead::read_reconstruction(const std::string& read_sequence_1d,
                                                                 uint32_t strand,
                                                                 std::vector<fast5::Basecall_Event>& f5_events,
                                                                 size_t k)
{
    // reconstruct sequence from event table
    assert(f_p and f_p->is_open());
    std::vector<EventRangeForBase> out_event_map;

    // initialize - one entry per read kmer
    uint32_t n_read_kmers = read_sequence_1d.size() - k + 1;
    out_event_map.resize(n_read_kmers);

    size_t curr_k_idx = 0;
    size_t curr_event_idx = 0;

    // Albacore represents `move` as how much sequence to add starting from the middle of the event kmer
    assert(k == 5);
    assert(curr_event_idx == 0);

    // Initialize the reconstructed sequence with the first 3-mer of the first k-mer
#if DEBUG_RECONSTRUCTION
    fprintf(stderr, "processing fast5: %s\n", this->fast5_path.c_str());
#endif

    // We keep track of the number of mismatches between the read sequence
    // and the event table. We skip duplicate mismatches to the same kmer
    // by tracking the kmen in previous_event_kmer
    std::string previous_event_kmer = "";
    size_t distinct_mismatches = 0;

    // Initialize
    out_event_map[curr_k_idx].indices[strand].start = curr_event_idx;
    curr_event_idx = 1;

    while(curr_event_idx < f5_events.size()) {
        int move = f5_events[curr_event_idx].move;
        if(move > 0) {
            // End current range
            out_event_map[curr_k_idx].indices[strand].stop = curr_event_idx - 1;

            // start new range
            curr_k_idx += move;

            // parse error
            if(curr_k_idx >= n_read_kmers) {
                break;
            }

            assert(curr_k_idx < out_event_map.size());
            out_event_map[curr_k_idx].indices[strand].start = curr_event_idx;
        }

        std::string inferred_kmer = read_sequence_1d.substr(curr_k_idx, k);
        std::string event_kmer = array2str(f5_events[curr_event_idx].model_state);

        // Check if the k-mer in the read sequene matches that in the event table
        // If not we increment the counter to predict whether the parse was succesfull
        if(inferred_kmer != event_kmer) {
            distinct_mismatches += event_kmer != previous_event_kmer;
#if DEBUG_RECONSTRUCTION
            std::string long_context = read_sequence_1d.substr(curr_k_idx - k, 3*k);
            fprintf(stderr, "[reconstruction] k:%zu e:%zu %s %s match? %zu move: %zu context: %s\n",
                curr_k_idx, curr_event_idx,
                inferred_kmer.c_str(), event_kmer.c_str(),
                inferred_kmer == event_kmer, f5_events[curr_event_idx].move, long_context.c_str());
#endif
        }

        curr_event_idx += 1;
        previous_event_kmer = event_kmer;
    }

    const double MISMATCH_THRESHOLD = 0.05;
    double mismatch_rate = (float)distinct_mismatches / n_read_kmers;
    if(mismatch_rate > MISMATCH_THRESHOLD || curr_k_idx >= n_read_kmers) {
        // poor read, skip
        g_unparseable_reads += 1;
        out_event_map.clear();
    } else {
        // good read, continue
        out_event_map[curr_k_idx].indices[strand].stop = curr_event_idx - 1;
    }

#if DEBUG_RECONSTRUCTION
    std::string classification = !out_event_map.empty() ? "goodread" : "badread";
    fprintf(stderr, "[final] %s %zu out of %zu kmers mismatch (%.2lf) classification: %s\n", this->fast5_path.c_str(), distinct_mismatches, n_read_kmers, mismatch_rate, classification.c_str());
#endif
    return out_event_map;
}

std::vector<EventRangeForBase> SquiggleRead::build_event_map_1d(const std::string& read_sequence_1d,
                                                                uint32_t strand,
                                                                std::vector<fast5::Basecall_Event>& f5_events,
                                                                size_t k)
{
    return read_reconstruction(read_sequence_1d, strand, f5_events, k);
}

void SquiggleRead::build_event_map_2d_r9()
{
    assert(f_p and f_p->is_open());
    assert(not basecall_group.empty());
    //
    // Build the map from read k-mers to events
    //
    auto event_alignments = f_p->get_basecall_alignment(basecall_group);
    assert(!read_sequence.empty());

    // R9 change: use k from the event table as this might not match the pore model
    uint32_t k = strnlen(event_alignments[0].kmer.data(), event_alignments[0].kmer.size());

    uint32_t n_read_kmers = read_sequence.size() - k + 1;
    base_to_event_map.resize(n_read_kmers);

    uint32_t read_kidx = 0;

    // The alignment format in the fast5 file is slightly bizarre in that it is
    // (template_idx, complement_idx, kmer) tuples. Some read kmers may not have a
    // tuple and some might have multiple tuples. We need to use the read kmer
    // sequences to work out which read base each entry is referring to
    uint32_t start_ea_idx = 0;
    uint32_t end_ea_idx = 0;
    //printf("Starting event map construction for read %s\n", read_name.c_str());
    while(start_ea_idx < event_alignments.size()) {

        uint32_t prev_kidx = read_kidx;

        // Advance the kmer index until we have found the read kmer
        // this tuple refers to
        while(read_kidx < n_read_kmers &&
              read_sequence.compare(read_kidx, k,
                                    array2str(event_alignments[start_ea_idx].kmer),
                                    0, k) != 0)
        {
            read_kidx += 1;
        }

        // Advance the event alignment end index to the last tuple
        // with the same kmer as the start of this range
        end_ea_idx = start_ea_idx;
        while(end_ea_idx < event_alignments.size() &&
                array2str(event_alignments[start_ea_idx].kmer).compare(0, k,
                          array2str(event_alignments[end_ea_idx].kmer), 0, k) == 0)
        {
            end_ea_idx += 1;
        }

        //printf("Base-to-event map kidx: %d %s event_tuple [%d %d]\n", read_kidx, read_sequence.substr(read_kidx, k).c_str(), start_ea_idx, end_ea_idx);
        EventRangeForBase& erfb =  base_to_event_map[read_kidx];
        for(uint32_t i = start_ea_idx; i < end_ea_idx; ++i) {

            auto const & eae = event_alignments[i];

            for(uint32_t si = 0; si <= 1; ++si) {
                int incoming_idx = si == 0 ? eae.template_index : eae.complement_index;

                // if a strand couldn't be loaded (typically because calibration failed) ignore
                // the events for the strand by setting to -1
                incoming_idx = events[si].empty() ? -1 : incoming_idx;

                // no event for this strand, nothing to update
                if(incoming_idx == -1) {
                    continue;
                }
                if(erfb.indices[si].start == -1) {
                    erfb.indices[si].start = incoming_idx;
                }
                erfb.indices[si].stop = incoming_idx;

                assert(erfb.indices[si].start < (int)events[si].size());
                assert(erfb.indices[si].stop < (int)events[si].size());
            }
        }
        //printf("\t[%d %d] [%d %d]\n", erfb.indices[0].start, erfb.indices[0].stop, erfb.indices[1].start, erfb.indices[1].stop);
        start_ea_idx = end_ea_idx;
    }
}

//as above but with a hack to get around a metrichor bug...
void SquiggleRead::build_event_map_2d_r7()
{
    assert(f_p and f_p->is_open());
    assert(not basecall_group.empty());
    //
    // Build the map from read k-mers to events
    //
    auto event_alignments = f_p->get_basecall_alignment(basecall_group);
    assert(!read_sequence.empty());

    const uint32_t k = this->get_model_k(T_IDX);
    assert(this->get_model_k(C_IDX) == k);
    assert(k == 5 || k == 6);

    uint32_t n_read_kmers = read_sequence.size() - k + 1;
    base_to_event_map.resize(n_read_kmers);

    uint32_t read_kidx = 0;

    // The alignment format in the fast5 file is slightly bizarre in that it is
    // (template_idx, complement_idx, kmer) tuples. Some read kmers may not have a
    // tuple and some might have multiple tuples. We need to use the read kmer
    // sequences to work out which read base each entry is referring to
    uint32_t start_ea_idx = 0;
    uint32_t end_ea_idx = 0;

    while(start_ea_idx < event_alignments.size()) {

hack:
        uint32_t prev_kidx = read_kidx;

        // Advance the kmer index until we have found the read kmer
        // this tuple refers to
        while(read_kidx < n_read_kmers &&
              read_sequence.compare(read_kidx, k,
                                    array2str(event_alignments[start_ea_idx].kmer),
                                    0, k) != 0)
        {
            read_kidx += 1;
        }
 
        // In the most recent version of metrichor occasionally
        // a kmer will be present in the alignment table
        // that is not in the 2D read. This awful hack
        // will skip such k-mers. It is not a long-term
        // solution, only until metrichor is fixed.
        if(read_kidx - prev_kidx > 10) {
            start_ea_idx += 1;
            read_kidx = prev_kidx;
            goto hack;
        }

        // Advance the event alignment end index to the last tuple
        // with the same kmer as the start of this range
        end_ea_idx = start_ea_idx;
        while(end_ea_idx < event_alignments.size() &&
                array2str(event_alignments[start_ea_idx].kmer).compare(0, k,
                          array2str(event_alignments[end_ea_idx].kmer), 0, k) == 0)
        {
            end_ea_idx += 1;
        }

        //printf("Base-to-event map kidx: %d %s event_tuple [%d %d]\n", read_kidx, read_sequence.substr(read_kidx, k).c_str(), start_ea_idx, end_ea_idx);
        EventRangeForBase& erfb =  base_to_event_map[read_kidx];
        for(uint32_t i = start_ea_idx; i < end_ea_idx; ++i) {

            auto const & eae = event_alignments[i];

            for(uint32_t si = 0; si <= 1; ++si) {
                uint32_t incoming_idx = si == 0 ? eae.template_index : eae.complement_index;

                // no event for this strand, nothing to update
                if(incoming_idx == -1)
                    continue;

                if(erfb.indices[si].start == -1) {
                    erfb.indices[si].start = incoming_idx;
                }
                erfb.indices[si].stop = incoming_idx;

                assert(erfb.indices[si].start < events[si].size());
                assert(erfb.indices[si].stop < events[si].size());
            }
        }
        //printf("\t[%d %d] [%d %d]\n", erfb.indices[0].start, erfb.indices[0].stop, erfb.indices[1].start, erfb.indices[1].stop);
        start_ea_idx = end_ea_idx;
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
    return sample_time - sample_start_time;
}

//
std::vector<float> SquiggleRead::get_scaled_samples_for_event(size_t strand_idx, size_t event_idx) const
{
    double event_start_time = this->events[strand_idx][event_idx].start_time;
    double event_duration = this->events[strand_idx][event_idx].duration;

    size_t start_idx = this->get_sample_index_at_time(event_start_time * this->sample_rate);
    size_t end_idx = this->get_sample_index_at_time((event_start_time + event_duration) * this->sample_rate);

    std::vector<float> out;
    for(size_t i = start_idx; i < end_idx; ++i) {
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

void SquiggleRead::detect_pore_type()
{
    assert(f_p and f_p->is_open());
    if (f_p->have_basecall_model(0))
    {
        pore_type = PT_R7;
    }
    else
    {
        pore_type = PT_R9;
    }
}

bool SquiggleRead::check_basecall_group() const
{
    assert(f_p and f_p->is_open());
    assert(not basecall_group.empty());
    if (not f_p->have_basecall_seq(read_type, basecall_group)) return false;
    if ((read_type == SRT_2D or read_type == SRT_TEMPLATE)
        and not f_p->have_basecall_events(0, basecall_group)) return false;
    if ((read_type == SRT_2D or read_type == SRT_COMPLEMENT)
        and not f_p->have_basecall_events(1, basecall_group)) return false;
    return true;
}

//
bool SquiggleRead::is_extract_read_name(std::string& name) const
{
    // albacore read names are uuids with hex characters separated
    // by underscores. If the read name contains three colon-delimited fields
    // we infer it is output by nanopolish extract. 
    return std::count(name.begin(), name.end(), ':') == 2;
}

void SquiggleRead::detect_basecall_group()
{
    assert(f_p and f_p->is_open());
    basecall_group.clear();
    // option 1: parse read name as written by extract:
    //   <uuid>:<basecall_group>:<template|complement|2d>
    do
    {
        auto i = read_name.find_first_of(':');
        if (i == std::string::npos) break;
        // 0..i : uuid
        ++i;
        auto j = read_name.find_first_of(':', i);
        if (j == std::string::npos) break;
        // i..j : basecall group
        basecall_group = read_name.substr(i, j - i);
        i = j + 1;
        j = read_name.find_first_of(':', i);
        if (j == std::string::npos) j = read_name.size();
        // i..j : read type
        auto rt = read_name.substr(i, j - i);
        if (not (rt == "2d"
                 or rt == "template"
                 or rt == "complement")) break;
        read_type = (rt == "2d" ? SRT_2D : (rt == "template"? SRT_TEMPLATE : SRT_COMPLEMENT));
        // check detection results
        if (not check_basecall_group())
        {
            std::cerr << "SquiggleRead: basecall group detection failed for name written by extract\n"
                      << "SquiggleRead: file name [" << fast5_path << "]\n"
                      << "SquiggleRead: read name [" << read_name << "]\n"
                      << "SquiggleRead: basecall group [" << basecall_group << "]\n"
                      << "SquiggleRead: read_type [" << int(read_type) << "]\n"
                      << "SquiggleRead: please submit a bug report about this\n";
            exit(EXIT_FAILURE);
        }
        return;
    }
    while (false);
    // option 2: find a Fastq record matching the read name
    auto gr_l = f_p->get_basecall_group_list();
    for (const auto& g : gr_l)
    {
        for (unsigned st = 0; st < 3; ++st)
        {
            if (not f_p->have_basecall_fastq(st, g)) continue;
            auto fq = f_p->get_basecall_fastq(st, g);
            auto fq_a = f_p->split_fq(fq);
            auto p = std::mismatch(read_name.begin(), read_name.end(), fq_a[0].begin());
            if (p.first != read_name.end()) continue;
            if (not basecall_group.empty())
            {
                std::cerr << "SquiggleRead: basecall group detection failed:\n"
                          << "SquiggleRead: multiple Fastq records match read name\n"
                          << "SquiggleRead: file name [" << fast5_path << "]\n"
                          << "SquiggleRead: read name [" << read_name << "]\n"
                          << "SquiggleRead: giving up; re-extract reads with 'nanopolish extract'\n";
                exit(EXIT_FAILURE);
            }
            basecall_group = g;
            read_type = (st == 2? SRT_2D : (st == 0? SRT_TEMPLATE : SRT_COMPLEMENT));
        }

        if(not basecall_group.empty()) {
            break;
        }
    }
    if (basecall_group.empty() or not check_basecall_group())
    {
        std::cerr << "SquiggleRead: basecall group detection failed:\n"
                  << "SquiggleRead: file name [" << fast5_path << "]\n"
                  << "SquiggleRead: read name [" << read_name << "]\n"
                  << "SquiggleRead: basecall group [" << basecall_group << "]\n"
                  << "SquiggleRead: read_type [" << (not basecall_group.empty()? int(read_type) : -1) << "]\n"
                  << "SquiggleRead: giving up; re-extract reads with 'nanopolish extract'\n";
        exit(EXIT_FAILURE);
    }
}
