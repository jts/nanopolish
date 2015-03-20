//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_squiggle_read -- Class holding a squiggle (event)
// space nanopore read
//
#include "nanopolish_squiggle_read.h"
#include "fast5/src/fast5.hpp"

//
SquiggleRead::SquiggleRead(const std::string& name, const std::string& path) :
    read_name(name)
{
    load_from_fast5(path);
}

//
void SquiggleRead::load_from_fast5(const std::string& fast5_path)
{
    printf("Loading %s\n", fast5_path.c_str());

    fast5::File* f_p;
    f_p = new fast5::File(fast5_path);
    assert(f_p->is_open());

    //
    std::cout << "file_version=" << f_p->file_version() << std::endl;
    std::cout << "basecall_version=" << f_p->basecall_version() << std::endl;
    std::cout << "eventdetection_version=" << f_p->eventdetection_version() << std::endl;
    std::cout << "sequences_version=" << f_p->sequences_version() << std::endl;

    // Load PoreModel for both strands
    for (size_t si = 0; si < 2; ++si) {

        std::vector<fast5::Model_Entry> model = f_p->get_model(si);
        assert(model.size() == 1024);
        assert(strcmp(model[0].kmer, "AAAAA") == 0);
        assert(strcmp(model[1023].kmer, "TTTTT") == 0);

        // Copy into the pore model for this read
        for(size_t mi = 0; mi < model.size(); ++mi) {
            const fast5::Model_Entry& curr = model[mi];
            pore_model[si].state[mi] = { curr.level_mean, curr.level_stdv, curr.sd_mean, curr.sd_stdv };
        }

        // Load the scaling parameters for the pore model
        fast5::Model_Parameters params = f_p->get_model_parameters(si);
        pore_model[si].drift = params.drift;
        pore_model[si].scale = params.scale;
        pore_model[si].scale_sd = params.scale_sd;
        pore_model[si].shift = params.shift;
        pore_model[si].var = params.var;
        pore_model[si].var_sd = params.var_sd;

    }
    
    printf("Template AAAAA: %.2lf %.2lf\n",  pore_model[0].state[0].level_mean, pore_model[0].state[0].level_stdv);
    printf("Complemt AAAAA: %.2lf %.2lf\n",  pore_model[1].state[0].level_mean, pore_model[1].state[0].level_stdv);

    // Load events for both strands
    for (size_t si = 0; si < 2; ++si) {
        std::vector<fast5::Event_Entry> f5_events = f_p->get_events(si);
        
        // copy events
        events[si].resize(f5_events.size());
        for(size_t ei = 0; ei < f5_events.size(); ++ei) {
            const fast5::Event_Entry& f5_event = f5_events[ei];
            events[si][ei] = { f5_event.mean, f5_event.stdv, f5_event.start, f5_event.length }; 
        }
    }

    printf("Loaded %zu template and %zu complement events\n", events[0].size(), events[1].size());
    delete f_p;
}
