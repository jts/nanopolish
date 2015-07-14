//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alignment_db -- abstraction for working
// with sets of reads/events aligned to a reference genome
//
#include <assert.h>
#include "nanopolish_alignment_db.h"
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"

// Various file handle and structures
// needed to traverse a bam file
struct BamHandles
{
    htsFile* bam_fh;
    bam1_t* bam_record;
    hts_itr_t* itr;
};

AlignmentDB::AlignmentDB(const std::string& reads_file,
                         const std::string& reference_file,
                         const std::string& sequence_bam,
                         const std::string& event_bam) :
                            m_fast5_name_map(reads_file),
                            m_reference_file(reference_file),
                            m_sequence_bam(sequence_bam),
                            m_event_bam(event_bam)
{
    _clear_region();
}

AlignmentDB::~AlignmentDB()
{
    _clear_region();
}
        
void AlignmentDB::load_region(const std::string& contig,
                              int start_position,
                              int stop_position)
{

    // load reference fai file
    faidx_t *fai = fai_load(m_reference_file.c_str());

    // Adjust end position to make sure we don't go out-of-range
    m_region_contig = contig;
    m_region_start = start_position;
    m_region_end = std::min(stop_position, faidx_seq_len(fai, contig.c_str()));

    // load the reference sequence for this region
    int fetched_len = 0;
    char* ref_segment = faidx_fetch_seq(fai, m_region_contig.c_str(), m_region_start, m_region_end, &fetched_len);
    m_region_ref_sequence = ref_segment;
    
    // load base-space alignments
    _load_sequence_by_region();

    // load event-space alignments
    _load_events_by_region();

    free(ref_segment);
    fai_destroy(fai);
}

void AlignmentDB::_clear_region()
{
    // Delete the SquiggleReads
    for(SquiggleReadMap::iterator iter = m_squiggle_read_map.begin();
        iter != m_squiggle_read_map.end(); ++iter) 
    {
        delete iter->second;
    }
    m_squiggle_read_map.clear();
    m_sequence_records.clear();
    m_event_records.clear();

    m_region_contig = "";
    m_region_start = -1;
    m_region_end = -1;
}

BamHandles _initialize_bam_itr(const std::string& bam_filename,
                               const std::string& contig,
                               int start_position,
                               int stop_position)
{
    BamHandles handles;

    // load bam file
    handles.bam_fh = sam_open(bam_filename.c_str(), "r");
    assert(handles.bam_fh != NULL);

    // load bam index file
    std::string index_filename = bam_filename + ".bai";
    hts_idx_t* bam_idx = bam_index_load(index_filename.c_str());
    assert(bam_idx != NULL);

    // read the bam header to get the contig ID
    bam_hdr_t* hdr = sam_hdr_read(handles.bam_fh);
    int contig_id = bam_name2id(hdr, contig.c_str());
    
    // Initialize iteration
    handles.bam_record = bam_init1();
    handles.itr = sam_itr_queryi(bam_idx, contig_id, start_position, stop_position);

    hts_idx_destroy(bam_idx);
    bam_hdr_destroy(hdr);
    return handles;
}

void AlignmentDB::_load_sequence_by_region()
{
    assert(!m_region_contig.empty());
    assert(m_region_start >= 0);
    assert(m_region_end >= 0);

    BamHandles handles = _initialize_bam_itr(m_sequence_bam, m_region_contig, m_region_start, m_region_end);

    int result;
    while((result = sam_itr_next(handles.bam_fh, handles.itr, handles.bam_record)) >= 0) {
        SequenceAlignmentRecord seq_record;

        // copy sequence out of the record
        uint8_t* pseq = bam_get_seq(handles.bam_record);
        seq_record.sequence.resize(handles.bam_record->core.l_qseq);
        for(int i = 0; i < handles.bam_record->core.l_qseq; ++i) {
            seq_record.sequence[i] = seq_nt16_str[bam_seqi(pseq, i)];
        }
        
        // copy read base-to-reference alignment
        seq_record.aligned_bases = get_aligned_pairs(handles.bam_record);
        m_sequence_records.push_back(seq_record);
        
        printf("sequence_record[%zu] prefix: %s alignstart: [%d %d]\n", 
            m_sequence_records.size() - 1,
            m_sequence_records.back().sequence.substr(0, 20).c_str(),
            m_sequence_records.back().aligned_bases.front().ref_pos,
            m_sequence_records.back().aligned_bases.front().read_pos);
    }

    // cleanup
    sam_itr_destroy(handles.itr);
    bam_destroy1(handles.bam_record);
    sam_close(handles.bam_fh);
}

void AlignmentDB::_load_events_by_region()
{
    assert(!m_region_contig.empty());
    assert(m_region_start >= 0);
    assert(m_region_end >= 0);

    BamHandles handles = _initialize_bam_itr(m_event_bam, m_region_contig, m_region_start, m_region_end);

    int result;
    while((result = sam_itr_next(handles.bam_fh, handles.itr, handles.bam_record)) >= 0) {
        EventAlignmentRecord event_record;

        std::string full_name = bam_get_qname(handles.bam_record);
        
        // Check for the template/complement suffix
        size_t suffix_pos = 0;
        suffix_pos = full_name.find(".template");
        if(suffix_pos == std::string::npos) {
            suffix_pos = full_name.find(".complement");
        }
        assert(suffix_pos != std::string::npos);
        std::string read_name = full_name.substr(0, suffix_pos);
        std::string fast5_path = m_fast5_name_map.get_path(read_name);

        // Do we need to load this fast5 file?
        if(m_squiggle_read_map.find(read_name) == m_squiggle_read_map.end()) {
            m_squiggle_read_map[read_name] = new SquiggleRead(read_name, fast5_path);
        }
        event_record.sr = m_squiggle_read_map[read_name];

        // extract the event stride tag which tells us whether the
        // event indices are increasing or decreasing
        assert(bam_aux_get(handles.bam_record, "ES") != NULL);
        int event_stride = bam_aux2i(bam_aux_get(handles.bam_record, "ES"));

        // copy event alignments
        event_record.aligned_events = get_aligned_pairs(handles.bam_record, event_stride);
        m_event_records.push_back(event_record);
        
        printf("event_record[%zu] name: %s stride: %d align bounds [%d %d] [%d %d]\n", 
            m_event_records.size() - 1,
            bam_get_qname(handles.bam_record),
            event_stride,
            m_event_records.back().aligned_events.front().ref_pos,
            m_event_records.back().aligned_events.front().read_pos,
            m_event_records.back().aligned_events.back().ref_pos,
            m_event_records.back().aligned_events.back().read_pos);
    }

    // cleanup
    sam_itr_destroy(handles.itr);
    bam_destroy1(handles.bam_record);
    sam_close(handles.bam_fh);
}

