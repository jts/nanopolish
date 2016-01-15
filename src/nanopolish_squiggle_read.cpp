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
#include "src/fast5.hpp"

//
SquiggleRead::SquiggleRead(const std::string& name, const std::string& path) :
    read_name(name),
    drift_correction_performed(false)
{
    load_from_fast5(path);

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
void SquiggleRead::load_from_fast5(const std::string& fast5_path)
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
    
    f_p->set_basecalled_group_id(group_id);

    // Automatically detect the read type from the name
    // Set the read type based on the suffix and load the sequence
    assert(!read_name.empty());
    if(ends_with(read_name, "_template")) {
        read_type = SRT_TEMPLATE;
        read_sequence = f_p->basecalled_1D(read_type);
    } else if(ends_with(read_name, "_complement")) {
        read_type = SRT_COMPLEMENT;
        read_sequence = f_p->basecalled_1D(read_type);
    } else {
        read_type = SRT_2D;
        read_sequence = f_p->basecalled_2D();
    }

    // Load PoreModel for both strands
    for (size_t si = 0; si < 2; ++si) {
        
        // Do we want to load this strand?
        if(! (read_type == SRT_2D || read_type == si) ) {
            continue;
        }

        // Load the pore model for this strand
        pore_model[si] = PoreModel( f_p, si );

        // initialize transition parameters
        parameters[si].initialize(pore_model[si].name);
    
        // Load the events for this strand
        std::vector<fast5::Event_Entry> f5_events = f_p->get_events(si);
        
        // copy events
        events[si].resize(f5_events.size());
        for(size_t ei = 0; ei < f5_events.size(); ++ei) {
            const fast5::Event_Entry& f5_event = f5_events[ei];
            events[si][ei] = { static_cast<float>(f5_event.mean), 
                               static_cast<float>(f5_event.stdv), 
                               static_cast<float>(f5_event.start), 
                               static_cast<float>(f5_event.length),
                               static_cast<float>(log(f5_event.stdv)) };
        }

        // build the map from k-mers to events if this is a 1D read
        if(read_type != SRT_2D) {
            build_event_map_1d(f_p, si, f5_events);
        }
    }

    // Both strands need to be loaded before the 2D event map is built
    if(read_type == SRT_2D) {
        build_event_map_2d(f_p);
    }

    delete f_p;
}

void SquiggleRead::build_event_map_1d(fast5::File* f_p, uint32_t strand, std::vector<fast5::Event_Entry>& f5_events)
{
    const uint32_t k = pore_model[T_IDX].k;
    
    // initialize - one entry per read kmer
    uint32_t n_read_kmers = read_sequence.size() - k + 1;
    base_to_event_map.resize(n_read_kmers);

    // The range for the first k-mer always starts at event 0
    assert(f5_events[0].move == 0);
    base_to_event_map[0].indices[strand].start = 0;

    size_t curr_k_idx = 0;
    for(size_t ei = 0; ei < f5_events.size(); ++ei) {
        const fast5::Event_Entry& f5_event = f5_events[ei];
        
        // Does this event correspond to a different k-mer than the previous one?
        if(f5_event.move > 0) {
            assert(ei != 0);

            // end the range for the current k-mer
            base_to_event_map[curr_k_idx].indices[strand].stop = ei - 1;
            curr_k_idx += f5_event.move;
            
            // start the range for the next kmer
            base_to_event_map[curr_k_idx].indices[strand].start = ei; 
        }

        // ensure the traversal matches the read sequence
        assert(read_sequence.compare(curr_k_idx, k, f5_event.model_state) == 0);
    }
}

void SquiggleRead::build_event_map_2d(fast5::File* f_p)
{
    //
    // Build the map from read k-mers to events
    //
    
    std::vector<fast5::Event_Alignment_Entry> event_alignments = f_p->get_event_alignments();
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
              strncmp(event_alignments[start_ea_idx].kmer, 
                     read_sequence.c_str() + read_kidx, k) != 0) {
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
                strcmp(event_alignments[start_ea_idx].kmer,
                    event_alignments[end_ea_idx].kmer) == 0) {
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

void SquiggleRead::replace_models(const ModelMap& map) {
    
    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {
        std::string curr_model = this->pore_model[strand_idx].name;
        auto model_iter = map.find(curr_model);
        if(model_iter != map.end()) {
            this->pore_model[strand_idx].update_states( model_iter->second );
        } else {
            fprintf(stderr, "Error, alternative model %s not found\n", curr_model.c_str());
            exit(EXIT_FAILURE);
        }
    }
}
