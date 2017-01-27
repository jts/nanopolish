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
#include <fast5.hpp>

//#define DEBUG_MODEL_SELECTION 1

//
SquiggleRead::SquiggleRead(const std::string& name, const std::string& path, const uint32_t flags) :
    read_name(name),
    pore_type(PT_UNKNOWN),
    fast5_path(path),
    drift_correction_performed(false),
    f_p(nullptr)
{
    load_from_fast5(flags);

    // perform drift correction and other scalings
    transform();
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
    uint32_t k = pore_model[strand].k;

    int stop_before = std::max(0, k_idx - 1000);
    int stop_after = std::min(k_idx + 1000, (int32_t)read_sequence.size() - (int32_t)k + 1);

    int event_before = get_next_event(k_idx, stop_before, -1, strand);
    int event_after = get_next_event(k_idx, stop_after, 1, strand);

    // TODO: better selection of "best" event to return
    if(event_before == -1)
        return event_after;
    return event_before;
}

//
void SquiggleRead::transform()
{
    for (size_t si = 0; si < 2; ++si) {
        for(size_t ei = 0; ei < events[si].size(); ++ei) {

            SquiggleEvent& event = events[si][ei];

            // correct level by drift
            double time = event.start_time - events[si][0].start_time;
            event.mean -= (time * pore_model[si].drift);
        }
    }

    drift_correction_performed = true;
}

//
void SquiggleRead::load_from_fast5(const uint32_t flags)
{
    f_p = new fast5::File(fast5_path);
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
        std::vector<fast5::Event_Entry> f5_events = f_p->get_basecall_events(si, basecall_group);

        // copy events
        events[si].resize(f5_events.size());
        std::vector<double> p_model_states;

        for(size_t ei = 0; ei < f5_events.size(); ++ei) {
            const fast5::Event_Entry& f5_event = f5_events[ei];

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
                                                      f_p->get_basecall_group_1d(basecall_group));
        event_maps_1d[si] = build_event_map_1d(read_sequences_1d[si], si, f5_events);

        // run version-specific load
        if(pore_type == PT_R7) {
            _load_R7(si);
        } else {
            _load_R9(si, read_sequences_1d[si], event_maps_1d[si], p_model_states, flags);
        }
    }

    // Build the map from k-mers of the read sequence to events
    if(read_type == SRT_2D) {
        if(pore_model[0].metadata.is_r9()) {
            build_event_map_2d_r9();
        } else {
            assert(pore_model[0].metadata.is_r7());
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
        fast5::Channel_Id_Parameters channel_params = f_p->get_channel_id_params();
        sample_rate = channel_params.sampling_rate;
    }

    delete f_p;
    f_p = nullptr;
}

void SquiggleRead::_load_R7(uint32_t si)
{
    assert(f_p and f_p->is_open());
    assert(not basecall_group.empty());
    // Load the pore model for this strand
    pore_model[si] = PoreModel( f_p, si, basecall_group );

    // initialize transition parameters
    parameters[si].initialize(get_model_metadata_from_name(pore_model[si].name));
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
    std::vector<EventAlignment> alignment =
        get_eventalignment_for_1d_basecalls(read_sequence_1d, event_map_1d, 5, si);

    // JTS Hack: blacklist bad k-mer and filter out events with low p_model_state
    double keep_fraction = 0.75;
    std::vector<double> sorted_p_model_states = p_model_states;
    std::sort(sorted_p_model_states.begin(), sorted_p_model_states.end());
    double p_model_state_threshold = sorted_p_model_states[sorted_p_model_states.size() * (1 - keep_fraction)];

    std::string blacklist_kmer = "CCTAG";
    std::vector<EventAlignment> filtered;

    assert(p_model_states.size() == events[si].size());

    size_t total_kmers = 0;
    std::string prev_kmer = "";

    for(const auto& ea : alignment) {
        if((!ea.rc && ea.ref_kmer == blacklist_kmer) ||
           (ea.rc && ea.ref_kmer == gDNAAlphabet.reverse_complement(blacklist_kmer)))
        {
            continue;
        }

        if(ea.ref_kmer != prev_kmer) {
            prev_kmer = ea.ref_kmer;
            total_kmers++;
        }

        if(p_model_states[ea.event_idx] < p_model_state_threshold)
            continue;

        filtered.push_back(ea);
    }

    events_per_base[si] = ((float)alignment.size() / total_kmers);

    // Load the pore model (if requested) and calibrate it
    if( (flags & SRF_NO_MODEL) == 0) {

        auto config = f_p->get_basecall_config(basecall_group);
        auto mt_iter = config.find("general/model_type");
        std::string kit = "r9_250bps"; //
        if(mt_iter != config.end()) {
            std::string mt = mt_iter->second;

            // all 250bps data should use this model (according to ONT see
            // https://github.com/nanoporetech/kmer_models/issues/3)
            if(mt == "r9_250bps_nn" || mt == "r9_250bps" || mt == "r94_250bps") {
                // intentially blank, default as above
            } else if(mt == "r94_450bps" || mt == "r9_450bps") {
                kit = "r9.4_450bps";
            } else {
                fprintf(stderr, "Unknown model type string: %s, please report on github.\n", mt.c_str());
                exit(1);
            }
        } else {
            WARN_ONCE("Warning: model type could not be detected from fast5, defaulting to R9-250bps")
        }

        std::string alphabet = "nucleotide"; // always calibrate with the nucleotide alphabet

        // For the template strad we only have one candidate model
        // For complement we need to select between the two possible models
        std::vector<const PoreModel*> candidate_models;
        if(si == 0) {
            candidate_models.push_back(&PoreModelSet::get_model(kit, alphabet, "template", calibration_k));
        } else {
            for(const std::string& cmn : { "complement.pop1", "complement.pop2" } ) {
                if(PoreModelSet::has_model(kit, alphabet, cmn, calibration_k)) {
                   candidate_models.push_back(&PoreModelSet::get_model(kit, alphabet, cmn, calibration_k));
                }
            }
        }

        PoreModel best_model;
        double best_model_var = INFINITY;

        for(size_t model_idx = 0; model_idx < candidate_models.size(); model_idx++) {

            pore_model[si] = *candidate_models[model_idx];

            // Initialize to default scaling parameters
            pore_model[si].shift = 0.0;
            pore_model[si].scale = 1.0;
            pore_model[si].drift = 0.0;
            pore_model[si].var = 1.0;
            pore_model[si].scale_sd = 1.0;
            pore_model[si].var_sd = 1.0;

            pore_model[si].bake_gaussian_parameters();

            // run recalibration to get the best set of scaling parameters and the residual
            // between the (scaled) event levels and the model
            bool calibrated = recalibrate_model(*this, si, filtered, pore_model[si].pmalphabet, true, false);
            if(calibrated) {
                if(pore_model[si].var < best_model_var) {
                    best_model_var = pore_model[si].var;
                    best_model = pore_model[si];
                }
            }

#ifdef DEBUG_MODEL_SELECTION
            fprintf(stderr, "[calibration] read: %s strand: %zu model_idx: %zu "
                             "scale: %.2lf shift: %.2lf drift: %.5lf var: %.2lf\n",
                                    read_name.substr(0, 6).c_str(), si, model_idx, pore_model[si].scale,
                                    pore_model[si].shift, pore_model[si].drift, pore_model[si].var);
#endif
        }

        if(best_model_var != INFINITY) {
#ifdef DEBUG_MODEL_SELECTION
            fprintf(stderr, "[calibration] selected model with var %.4lf\n", best_model_var);
#endif
            pore_model[si] = best_model;

            // Save calibration parameters
            double shift = pore_model[si].shift;
            double scale = pore_model[si].scale;
            double drift =  pore_model[si].drift;
            double var = pore_model[si].var;
            double scale_sd = pore_model[si].scale_sd;
            double var_sd = pore_model[si].var_sd;

            // Replace model
            PoreModel final_model = PoreModelSet::get_model(kit,
                                                            alphabet,
                                                            best_model.metadata.get_strand_model_name(),
                                                            final_model_k);
            pore_model[si] = final_model;

            // Copy calibration params
            pore_model[si].shift = shift;
            pore_model[si].scale = scale;
            pore_model[si].drift = drift;
            pore_model[si].var = var;
            pore_model[si].scale_sd = scale_sd;
            pore_model[si].var_sd = var_sd;

            // Initialize gaussian params
            pore_model[si].bake_gaussian_parameters();

            // initialize transition parameters
            parameters[si].initialize(pore_model[si].metadata);
        } else {
            // could not find a model for this strand, discard it
            events[si].clear();
        }
    }
}

std::vector<EventRangeForBase> SquiggleRead::build_event_map_1d(const std::string& read_sequence_1d,
                                                                uint32_t strand,
                                                                std::vector<fast5::Event_Entry>& f5_events)
{
    assert(f_p and f_p->is_open());
    std::vector<EventRangeForBase> out_event_map;
    const uint32_t k = pore_model[strand].k;
    assert(f5_events.size() == events[strand].size());

    // initialize - one entry per read kmer
    uint32_t n_read_kmers = read_sequence_1d.size() - k + 1;
    out_event_map.resize(n_read_kmers);

    // The range for the first k-mer always starts at event 0
    assert(f5_events[0].move == 0);
    out_event_map[0].indices[strand].start = 0;

    size_t curr_k_idx = 0;
    for(size_t ei = 0; ei < f5_events.size(); ++ei) {
        const fast5::Event_Entry& f5_event = f5_events[ei];

        // Does this event correspond to a different k-mer than the previous one?
        if(f5_event.move > 0) {
            assert(ei != 0);

            // end the range for the current k-mer
            out_event_map[curr_k_idx].indices[strand].stop = ei - 1;
            curr_k_idx += f5_event.move;

            // start the range for the next kmer
            out_event_map[curr_k_idx].indices[strand].start = ei;
        }

        assert(read_sequence_1d.compare(curr_k_idx, k,
                                        array2str(f5_event.model_state), 0, k) == 0);
    }

    // end the last range
    out_event_map[curr_k_idx].indices[strand].stop = events[strand].size() - 1;
    assert(out_event_map[curr_k_idx].indices[strand].start <= out_event_map[curr_k_idx].indices[strand].stop);
    return out_event_map;
}

void SquiggleRead::build_event_map_2d_r9()
{
    assert(f_p and f_p->is_open());
    assert(not basecall_group.empty());
    //
    // Build the map from read k-mers to events
    //
    std::vector<fast5::Event_Alignment_Entry> event_alignments = f_p->get_basecall_event_alignment(basecall_group);
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

            fast5::Event_Alignment_Entry& eae = event_alignments[i];

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
    std::vector<fast5::Event_Alignment_Entry> event_alignments = f_p->get_basecall_event_alignment(basecall_group);
    assert(!read_sequence.empty());

    const uint32_t k = pore_model[T_IDX].k;
    assert(pore_model[C_IDX].k == k);

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

            fast5::Event_Alignment_Entry& eae = event_alignments[i];

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

void SquiggleRead::replace_models(const std::string& kit_name, const std::string& alphabet, size_t k)
{
    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {

        // only replace this model if the strand was loaded
        if( !(read_type == SRT_2D || read_type == strand_idx) || !has_events_for_strand(strand_idx)) {
            continue;
        }

        PoreModel incoming_model =
            PoreModelSet::get_model(kit_name, 
                                    alphabet, 
                                    this->pore_model[strand_idx].metadata.get_strand_model_name(),
                                    k);
        replace_model(strand_idx, incoming_model);
    }
}

void SquiggleRead::replace_model(size_t strand_idx, const PoreModel& model)
{
    WARN_ONCE("Copy all model data in replace");
    this->pore_model[strand_idx].metadata = model.metadata;
    this->pore_model[strand_idx].k = model.k;
    this->pore_model[strand_idx].pmalphabet = model.pmalphabet;
    this->pore_model[strand_idx].update_states( model );
}

void SquiggleRead::replace_model(size_t strand_idx, const std::string& model_type)
{
    assert(false);
#if 0
    PoreModel incoming_model =
        PoreModelSet::get_model(model_type, this->pore_model[strand_idx].metadata.get_short_name());
    replace_model(strand_idx, incoming_model);
#endif
}

// Return a vector of eventalignments for the events that made up the 2D basecalls in the read
std::vector<EventAlignment> SquiggleRead::get_eventalignment_for_1d_basecalls(const std::string& read_sequence_1d,
                                                                              const std::vector<EventRangeForBase>& base_to_event_map_1d,
                                                                              const size_t k,
                                                                              const size_t strand_idx) const
{
    std::vector<EventAlignment> alignment;

    const Alphabet* alphabet = this->pore_model[strand_idx].pmalphabet;
    size_t n_kmers = read_sequence_1d.size() - k + 1;
    size_t prev_kmer_rank = -1;

    for(size_t ki = 0; ki < n_kmers; ++ki) {
        IndexPair event_range_for_kmer = base_to_event_map_1d[ki].indices[strand_idx];

        // skip kmers without events
        if(event_range_for_kmer.start == -1)
            continue;

        for(size_t event_idx = event_range_for_kmer.start;
            event_idx <= event_range_for_kmer.stop; event_idx++)
        {
            assert(event_idx < this->events[strand_idx].size());

            // since we use the 1D read seqence here we never have to reverse complement
            std::string kmer = read_sequence_1d.substr(ki, k);
            size_t kmer_rank = alphabet->kmer_rank(kmer.c_str(), k);

            EventAlignment ea;
            // ref data
            ea.ref_name = "read"; // not needed
            ea.ref_kmer = kmer;
            ea.ref_position = ki;
            ea.read_idx = -1; // not needed
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
        double scaled_s = s - this->pore_model[strand_idx].shift;
        assert(curr_sample_time >= (this->sample_start_time / this->sample_rate));
        scaled_s -= (curr_sample_time - (this->sample_start_time / this->sample_rate)) * this->pore_model[strand_idx].drift;
        scaled_s /= this->pore_model[strand_idx].scale;
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
