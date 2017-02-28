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
#include "nanopolish_variant.h"

#define MAX_EVENT_TO_BP_RATIO 20

// structs
struct SequenceAlignmentRecord
{
    SequenceAlignmentRecord(const bam1_t* record);

    std::string read_name;
    std::string sequence;
    std::vector<AlignedPair> aligned_bases;
    uint8_t rc; // with respect to reference genome
};

struct EventAlignmentRecord
{
    EventAlignmentRecord() {}
    EventAlignmentRecord(SquiggleRead* sr,
                         const int strand_idx,
                         const SequenceAlignmentRecord& seq_record);

    SquiggleRead* sr;
    uint8_t rc; // with respect to reference genome
    uint8_t strand; // 0 = template, 1 = complement
    int8_t stride; // whether event indices increase or decrease along the reference
    std::vector<AlignedPair> aligned_events;
};

// typedefs
typedef std::map<std::string, SquiggleRead*> SquiggleReadMap;

class AlignmentDB
{
    public:
        AlignmentDB(const std::string& reads_file,
                    const std::string& reference_file,
                    const std::string& sequence_bam,
                    const std::string& event_bam,
                    const bool calibrate_reads = false);

        ~AlignmentDB();

        void load_region(const std::string& contig,
                         int start_position,
                         int stop_position);
        
        const std::string& get_reference() const { return m_region_ref_sequence; }

        std::string get_reference_substring(const std::string& contig,
                                            int start_position,
                                            int stop_position) const;

        std::vector<std::string> get_read_substrings(const std::string& contig,
                                                     int start_position,
                                                     int stop_position) const;

        std::vector<HMMInputData> get_event_subsequences(const std::string& contig,
                                                         int start_position,
                                                         int stop_position) const;

        std::vector<HMMInputData> get_events_aligned_to(const std::string& contig, int position) const;

        std::vector<Variant> get_variants_in_region(const std::string& contig,
                                                    int start_position,
                                                    int stop_position,
                                                    double min_frequency,
                                                    int min_depth) const;

        const std::vector<EventAlignmentRecord>& get_eventalignment_records() const { return m_event_records; }

        // reference metadata
        std::string get_region_contig() const { return m_region_contig; }
        int get_region_start() const { return m_region_start; }
        int get_region_end() const { return m_region_end; }
        
        void set_alternative_model_type(const std::string model_type_string) { m_model_type_string = model_type_string; }
        
        // Search the vector of AlignedPairs using lower_bound/upper_bound
        // and the input reference coordinates. If the search succeeds,
        // set read_start/read_stop to be the read_pos of the bounding elements
        // and return true. 
        static bool _find_by_ref_bounds(const std::vector<AlignedPair>& pairs,
                                 int ref_start,
                                 int ref_stop,
                                 int& read_start,
                                 int& read_stop);

        static bool _find_iter_by_ref_bounds(const std::vector<AlignedPair>& pairs,
                                      int ref_start,
                                      int ref_stop,
                                      AlignedPairConstIter& start_iter,
                                      AlignedPairConstIter& stop_iter);
    private:
        
        void _load_sequence_by_region();
        void _load_events_by_region();
        void _load_events_by_region_from_bam();
        void _load_events_by_region_from_read();
        void _load_squiggle_read(const std::string& read_name);

        void _clear_region();

        void _debug_print_alignments();

        std::vector<EventAlignment> _build_event_alignment(const EventAlignmentRecord& event_record) const;


        //
        // data
        //
        std::string m_reference_file;
        std::string m_sequence_bam;
        std::string m_event_bam;

        // parameters
        bool m_calibrate_on_load;

        // loaded region
        std::string m_region_ref_sequence;
        std::string m_region_contig;
        int m_region_start;
        int m_region_end;

        // cached alignments for a region
        Fast5Map m_fast5_name_map;
        std::vector<SequenceAlignmentRecord> m_sequence_records;
        std::vector<EventAlignmentRecord> m_event_records;
        SquiggleReadMap m_squiggle_read_map;
        std::string m_model_type_string;
};

#endif
