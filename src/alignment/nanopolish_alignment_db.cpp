//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alignment_db -- abstraction for working
// with sets of reads/events aligned to a reference genome
//
#include "nanopolish_alignment_db.h"

AlignmentDB::AlignmentDB(const std::string& sequence_bam,
                         const std::string& event_bam) :
                            m_sequence_bam(sequence_bam),
                            m_event_bam(event_bam)
{

}
