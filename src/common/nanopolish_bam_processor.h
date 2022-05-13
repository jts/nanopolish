//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_bam_processor -- framework for iterating
// over a bam file and running an arbitrary function
// on each aligned read in parallel
//
#ifndef NANOPOLISH_BAM_PROCESSOR_H
#define NANOPOLISH_BAM_PROCESSOR_H

#include <functional>
#include <string>
#include <vector>
#include "htslib/hts.h"
#include "htslib/sam.h"

typedef std::function<void(const bam_hdr_t* hdr,
          const bam1_t* record,
          size_t read_idx,
          int region_start,
          int region_end)> bam_record_function;

typedef std::function<void(const bam_hdr_t* hdr,
          const std::vector<bam1_t*> batch)> batch_function;

class BamProcessor
{

    public:
        BamProcessor(const std::string& bam_filename, 
                     const std::string& region,
                     const int num_threads,
                     const int batch_size=512);

        ~BamProcessor();

        const bam_hdr_t* get_bam_header() const { return m_hdr; }

        // place a limit on the number of reads to process before stopping
        void set_max_reads(size_t max) { m_max_reads = max; }

        // place a limit on the minimum mapping quality
        void set_min_mapping_quality(size_t min_mapq) { m_min_mapping_quality = min_mapq; }

        // process each record in parallel, using the input function
        void parallel_run(bam_record_function record_func, batch_function batch_func = nullptr);

    private:
        std::string m_bam_file;
        std::string m_region;
    
        htsFile* m_bam_fh;
        hts_idx_t* m_bam_idx;
        bam_hdr_t* m_hdr;

        int m_batch_size = 512;
        int m_num_threads = 1;
        size_t m_max_reads = -1;
        int m_min_mapping_quality = 0;
};

#endif
