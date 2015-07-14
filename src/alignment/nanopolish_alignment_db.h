//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alignment_db -- abstraction for working
// with sets of reads/events aligned to a reference genome
//
#ifndef ALIGNMENT_DB
#define ALIGNMENT_DB

#include <string>

class AlignmentDB
{
    public:
        AlignmentDB(const std::string& sequence_bam,
                    const std::string& event_bam);

    private:

        std::string m_sequence_bam;
        std::string m_event_bam;

};

#endif
