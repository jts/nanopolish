//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_bam_processor -- framework for iterating
// over a bam file and running an arbitrary function
// on each aligned read in parallel
//
#include "nanopolish_bam_processor.h"
#include "nanopolish_common.h"
#include <assert.h>
#include <omp.h>
#include <vector>
#include <hdf5.h>

BamProcessor::BamProcessor(const std::string& bam_file,
                           const std::string& region,
                           const int num_threads) :
                            m_bam_file(bam_file),
                            m_region(region),
                            m_num_threads(num_threads)

{
    // load bam file
    m_bam_fh = sam_open(m_bam_file.c_str(), "r");
    assert(m_bam_fh != NULL);

    // load bam index file
    std::string index_filename = m_bam_file + ".bai";
    m_bam_idx = bam_index_load(index_filename.c_str());
    if(m_bam_idx == NULL) {
        bam_index_error_exit(m_bam_file);
    }

    // read the bam header
    m_hdr = sam_hdr_read(m_bam_fh);
    assert(m_hdr != NULL);
}

BamProcessor::~BamProcessor()
{
    bam_hdr_destroy(m_hdr);
    sam_close(m_bam_fh);
    hts_idx_destroy(m_bam_idx);
}

void BamProcessor::parallel_run( std::function<void(const bam_hdr_t* hdr, 
                                           const bam1_t* record,
                                           size_t read_idx,
                                           int region_start,
                                           int region_end)> func)
{
    assert(m_bam_fh != NULL);
    assert(m_bam_idx != NULL);
    assert(m_hdr != NULL);
    hts_itr_t* itr;

    // If processing a region of the genome, pass clipping coordinates to the work function
    int clip_start = -1;
    int clip_end = -1;

    if(m_region.empty()) {
        itr = sam_itr_queryi(m_bam_idx, HTS_IDX_START, 0, 0);
    } else {
        fprintf(stderr, "[bam process] iterating over region: %s\n", m_region.c_str());
        itr = sam_itr_querys(m_bam_idx, m_hdr, m_region.c_str());
        hts_parse_reg(m_region.c_str(), &clip_start, &clip_end);
    }

#ifndef H5_HAVE_THREADSAFE
    if(m_num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif
    
    // store number of threads so we can restore it after we're done
    int prev_num_threads = omp_get_num_threads();
    omp_set_num_threads(m_num_threads);

    // Initialize iteration
    std::vector<bam1_t*> records(m_batch_size, NULL);
    for(size_t i = 0; i < records.size(); ++i) {
        records[i] = bam_init1();
    }

    int result;
    size_t num_reads_realigned = 0;
    size_t num_records_buffered = 0;

    do {
        assert(num_records_buffered < records.size());

        // read a record into the next slot in the buffer
        result = sam_itr_next(m_bam_fh, itr, records[num_records_buffered]);
        num_records_buffered += result >= 0;

        // realign if we've hit the max buffer size or reached the end of file
        if(num_records_buffered == records.size() || result < 0 || (num_records_buffered + num_reads_realigned == m_max_reads)) {
            #pragma omp parallel for
            for(size_t i = 0; i < num_records_buffered; ++i) {
                bam1_t* record = records[i];
                size_t read_idx = num_reads_realigned + i;
                if( (record->core.flag & BAM_FUNMAP) == 0) {
                    func(m_hdr, record, read_idx, clip_start, clip_end);
                }
            }

            num_reads_realigned += num_records_buffered;
            num_records_buffered = 0;
        }
    } while(result >= 0 && num_reads_realigned < m_max_reads);

    assert(num_records_buffered == 0);

    // restore number of threads
    omp_set_num_threads(prev_num_threads);
 
    // cleanup   
    for(size_t i = 0; i < records.size(); ++i) {
        bam_destroy1(records[i]);
    }

    sam_itr_destroy(itr);
}
