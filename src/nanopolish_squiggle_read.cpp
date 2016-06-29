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
#include "src/fast5.hpp"

//
SquiggleRead::SquiggleRead(const std::string& name, const std::string& path, const uint32_t flags) :
    read_name(name),
    drift_correction_performed(false)
{
    load_from_fast5(path, flags);

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
    uint32_t k = pore_model[T_IDX].k;

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
void SquiggleRead::load_from_fast5(const std::string& fast5_path, const uint32_t flags)
{
    fast5::File* f_p;
    this->fast5_path = fast5_path;

    f_p = new fast5::File(fast5_path);
    assert(f_p->is_open());

    // Check if an alternative analysis group is present in the read name
    int group_id = -1;
    /*
    size_t bc_2d_pos = read_name.find("Basecall_2D");
    if(bc_2d_pos != std::string::npos) {
        int ret = sscanf(read_name.substr(bc_2d_pos).c_str(), "Basecall_2D_%03d_2d", &group_id);
    }
    */
    // default to 0 group
    if(group_id == -1) {
        group_id = 0;
    }

    //f_p->set_basecalled_group_id(group_id);

    auto available_groups = f_p->get_basecall_group_list();

    // precedence: 2D_NNN, RNN_1D_NNN, 1D_NNN
    std::string basecall_group = "1D_000";
    std::string event_group = "1D_000";

    read_type = SRT_TEMPLATE;

    for(auto g : available_groups) {
        if(g == "2D_000") {
            basecall_group = g;
            // for 2D reads we still take events from the 1D group
            read_type = SRT_2D;
        } else if(g == "RNN_1D_000" && g != "2D_000") {
            basecall_group = g;
            event_group = g;
            read_type = SRT_TEMPLATE;
        }
    }

    read_sequence = f_p->get_basecall_seq(basecall_group, read_type);

    // Load PoreModel for both strands
    std::vector<EventRangeForBase> event_maps_1d[NUM_STRANDS];
    std::string read_sequences_1d[NUM_STRANDS];

    for (size_t si = 0; si < 2; ++si) {

        // Do we want to load this strand?
        if(! (read_type == SRT_2D || read_type == si) ) {
            continue;
        }

        // Initialize to default scaling parameters
        pore_model[si].shift = 0.0;
        pore_model[si].scale = 1.0;
        pore_model[si].drift = 0.0;
        pore_model[si].var = 1.0;
        pore_model[si].scale_sd = 1.0;
        pore_model[si].var_sd = 1.0;

        // R9 change: load pore model from external file

        // this flag should only be used when running train-poremodel-from-basecalls
        // in this case there is no default model to load so we skip this step
        if( (flags & SRF_NO_MODEL) == 0) {
            pore_model[si] = PoreModelSet::get_model("base", "t.007");
            pore_model[si].bake_gaussian_parameters();

            // initialize transition parameters
            parameters[si].initialize(pore_model[si].metadata);
        }

        // Load the events for this strand
        std::vector<fast5::Event_Entry> f5_events = f_p->get_basecall_events(event_group, si);

        // copy events
        events[si].resize(f5_events.size());
        for(size_t ei = 0; ei < f5_events.size(); ++ei) {
            const fast5::Event_Entry& f5_event = f5_events[ei];
            events[si][ei] = { static_cast<float>(f5_event.mean),
                               static_cast<float>(f5_event.stdv),
                               f5_event.start,
                               static_cast<float>(f5_event.length),
                               static_cast<float>(log(f5_event.stdv)) };
        }

        // we need the 1D event map and sequence to calculate calibration parameters
        // these will be copied into the member fields later if this is a 1D read,
        // or discarded if this is a 2D read

        // NB we use event_group in this call rather than basecall_group as we want the 1D basecalls that match the events
        read_sequences_1d[si] = f_p->get_basecall_seq(event_group, si == 0 ? SRT_TEMPLATE : SRT_COMPLEMENT);
        event_maps_1d[si] = build_event_map_1d(f_p, read_sequences_1d[si], si, f5_events);
        std::vector<EventAlignment> alignment =
            get_eventalignment_for_1d_basecalls(read_sequences_1d[si], event_maps_1d[si], 5, si);

        bool calibrated = recalibrate_model(*this, si, alignment, pore_model[si].pmalphabet, true);
        if(!calibrated) {
            events[si].clear();
        }
    }

    // If we detected a 2D read we have to check that both strands loaded correctly
    // or else we downgrade it to a template-only read. If it was the template
    // strand that didn't load this is ok, the read will just be ignored later.
    if(read_type == SRT_2D) {
        if(events[0].empty() || events[1].empty()) {
            read_type = SRT_TEMPLATE;
        }
    }

    // Build the map from k-mers of the read sequence to events
    if(read_type == SRT_2D) {
        build_event_map_2d(f_p, basecall_group);
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
}

std::vector<EventRangeForBase> SquiggleRead::build_event_map_1d(fast5::File* f_p,
                                                                const std::string& read_sequence_1d, 
                                                                uint32_t strand, 
                                                                std::vector<fast5::Event_Entry>& f5_events)
{
    std::vector<EventRangeForBase> out_event_map;
    const uint32_t k = pore_model[T_IDX].k;
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

void SquiggleRead::build_event_map_2d(fast5::File* f_p, const std::string& basecall_group)
{
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

void SquiggleRead::replace_models(const std::string& model_type)
{

    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {

        // only replace this model if the strand was loaded
        if(! (read_type == SRT_2D || read_type == strand_idx) ) {
            continue;
        }

        PoreModel incoming_model =
            PoreModelSet::get_model(model_type, this->pore_model[strand_idx].metadata.get_short_name());
        replace_model(strand_idx, incoming_model);
    }
}

void SquiggleRead::replace_model(size_t strand_idx, const std::string& model_type)
{
    PoreModel incoming_model =
        PoreModelSet::get_model(model_type, this->pore_model[strand_idx].metadata.get_short_name());
    replace_model(strand_idx, incoming_model);
}

void SquiggleRead::replace_model(size_t strand_idx, const PoreModel& model)
{
    this->pore_model[strand_idx].update_states( model );
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

            std::string kmer = strand_idx == T_IDX ?
                                read_sequence_1d.substr(ki, k) :
                                alphabet->reverse_complement(read_sequence_1d.substr(ki, k));
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

