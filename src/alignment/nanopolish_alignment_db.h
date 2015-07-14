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
#include <vector>
#include <map>
#include "nanopolish_anchor.h"

// structs
struct SequenceAlignmentRecord
{
    std::string sequence;
    std::vector<AlignedPair> aligned_bases;
};

struct EventAlignmentRecord
{
    SquiggleRead* sr;
    std::vector<AlignedPair> aligned_events;
};

// typedefs
typedef std::map<std::string, SquiggleRead*> SquiggleReadMap;

class AlignmentDB
{
    public:
        AlignmentDB(const std::string& reference_file,
                    const std::string& sequence_bam,
                    const std::string& event_bam);

        void load_region(const std::string& contig,
                         int start_position,
                         int stop_position);

    private:
        
        void _load_sequence_by_region();
        void _load_events_by_region();
        void _clear_region();

        std::string m_reference_file;
        std::string m_sequence_bam;
        std::string m_event_bam;

        // loaded region
        std::string m_region_ref_sequence;
        std::string m_region_contig;
        int m_region_start;
        int m_region_end;

        // cached alignments for a region
        std::vector<SequenceAlignmentRecord> m_sequence_records;
        std::vector<EventAlignmentRecord> m_event_records;
        SquiggleReadMap m_squiggle_read_map;
};

#endif
