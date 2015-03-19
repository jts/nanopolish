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
SquiggleRead::SquiggleRead(const std::string& fast5_path)
{
    load_from_fast5(fast5_path);
}

//
void SquiggleRead::load_from_fast5(const std::string& fast5_path)
{
    printf("Loading %s\n", fast5_path.c_str());

    fast5::File* f_p;
    f_p = new fast5::File(fast5_path);
    assert(f_p->is_open());
    std::cout << "file_version=" << f_p->file_version() << std::endl;
    std::cout << "basecall_version=" << f_p->basecall_version() << std::endl;
    std::cout << "eventdetection_version=" << f_p->eventdetection_version() << std::endl;
    std::cout << "sequences_version=" << f_p->sequences_version() << std::endl;
    delete f_p;
}
