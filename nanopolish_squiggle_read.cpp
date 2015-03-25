//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_squiggle_read -- Class holding a squiggle (event)
// space nanopore read
//
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"
#include "fast5/src/fast5.hpp"

//
SquiggleRead::SquiggleRead(const std::string& name, const std::string& path) :
    read_name(name),
    drift_correction_performed(false)
{
    load_from_fast5(path);

    // perform drift correction
    transform();

    // TODO: refactor
    khmm_parameters_initialize(parameters[0]);
    khmm_parameters_initialize(parameters[1]);
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
    int stop_after = std::min(k_idx + 1000, (int32_t)read_sequence.size() - K + 1);
    
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
void SquiggleRead::load_from_fast5(const std::string& fast5_path)
{
    printf("Loading %s\n", fast5_path.c_str());

    fast5::File* f_p;
    f_p = new fast5::File(fast5_path);
    assert(f_p->is_open());

    // Load PoreModel for both strands
    for (size_t si = 0; si < 2; ++si) {

        std::vector<fast5::Model_Entry> model = f_p->get_model(si);
        assert(model.size() == PORE_MODEL_STATES);
        assert(strcmp(model[0].kmer, "AAAAA") == 0);
        assert(strcmp(model[PORE_MODEL_STATES - 1].kmer, "TTTTT") == 0);

        // Copy into the pore model for this read
        for(size_t mi = 0; mi < model.size(); ++mi) {
            const fast5::Model_Entry& curr = model[mi];
            pore_model[si].state[mi] = { static_cast<float>(curr.level_mean), 
                                         static_cast<float>(curr.level_stdv), 
                                         static_cast<float>(curr.sd_mean),
                                         static_cast<float>(curr.sd_stdv) };
        }

        // Load the scaling parameters for the pore model
        fast5::Model_Parameters params = f_p->get_model_parameters(si);
        pore_model[si].drift = params.drift;
        pore_model[si].scale = params.scale;
        pore_model[si].scale_sd = params.scale_sd;
        pore_model[si].shift = params.shift;
        pore_model[si].var = params.var;
        pore_model[si].var_sd = params.var_sd;

        // apply shift/scale transformation to the pore model states
        pore_model[si].bake_gaussian_parameters();
    }
    
    // Load events for both strands
    for (size_t si = 0; si < 2; ++si) {
        std::vector<fast5::Event_Entry> f5_events = f_p->get_events(si);
        
        // copy events
        events[si].resize(f5_events.size());
        for(size_t ei = 0; ei < f5_events.size(); ++ei) {
            const fast5::Event_Entry& f5_event = f5_events[ei];
            events[si][ei] = { static_cast<float>(f5_event.mean), 
                               static_cast<float>(f5_event.stdv), 
                               static_cast<float>(f5_event.start), 
                               static_cast<float>(f5_event.length) }; 
        }
    }

    printf("Loaded %zu template and %zu complement events\n", events[0].size(), events[1].size());

    //
    // Load basecalled sequence
    //
    read_sequence = f_p->basecalled_2D();
    
    //
    // Build the map from read k-mers to events
    //
    std::vector<fast5::Event_Alignment_Entry> event_alignments = f_p->get_event_alignments();
    assert(!read_sequence.empty());

    uint32_t n_read_kmers = read_sequence.size() - K + 1;
    base_to_event_map.resize(n_read_kmers);

    uint32_t read_kidx = 0;

    // The alignment format in the fast5 file is slightly bizarre in that it is
    // (template_idx, complement_idx, kmer) tuples. Some read kmers may not have a
    // tuple and some might have multiple tuples. We need to use the read kmer
    // sequences to work out which read base each entry is referring to
    uint32_t start_ea_idx = 0;
    uint32_t end_ea_idx = 0;

    while(start_ea_idx < event_alignments.size()) {
        
        // Advance the kmer index until we have found the read kmer
        // this tuple refers to
        while(read_kidx < n_read_kmers && 
              strncmp(event_alignments[start_ea_idx].kmer, 
                     read_sequence.c_str() + read_kidx, K) != 0) {
            read_kidx += 1;
        }

        // Advance the event alignment end index to the last tuple
        // with the same kmer as the start of this range
        end_ea_idx = start_ea_idx;
        while(end_ea_idx < event_alignments.size() &&
              strcmp(event_alignments[start_ea_idx].kmer, 
                     event_alignments[end_ea_idx].kmer) == 0) {
            end_ea_idx += 1;
        }

        //printf("Base-to-event map kidx: %d %s event_tuple [%d %d]\n", read_kidx, read_sequence.substr(read_kidx, K).c_str(), start_ea_idx, end_ea_idx);
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
            }
        }
        //printf("\t[%d %d] [%d %d]\n", erfb.indices[0].start, erfb.indices[0].stop, erfb.indices[1].start, erfb.indices[1].stop);
        start_ea_idx = end_ea_idx;
    }

    delete f_p;
}
