//---------------------------------------------------------
// Copyright 2021 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_basemods -- utility functions for working with
// base modification data
//
#ifndef NANOPOLISH_BASEMODS_H
#define NANOPOLISH_BASEMODS_H

#include <stdio.h>
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "nanopolish_squiggle_read.h"

//
// Structs
//
struct OutputHandles
{
    FILE* site_writer = NULL;
    htsFile* bam_writer = NULL;

    void write_site_header() {
        // Write header
        fprintf(site_writer, "chromosome\tstrand\tstart\tend\tread_name\t"
                             "log_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\t"
                             "num_calling_strands\tnum_motifs\tsequence\n");
    }

    void close() {
        if(site_writer != stdout) {
            fclose(site_writer);
            site_writer = NULL;
        }

        if(bam_writer != NULL) {
            hts_close(bam_writer);
            bam_writer = NULL;
        }
    }
};

struct ScoredSite
{
    ScoredSite()
    {
        ll_unmethylated[0] = 0;
        ll_unmethylated[1] = 0;
        ll_methylated[0] = 0;
        ll_methylated[1] = 0;
        strands_scored = 0;
    }

    std::string chromosome;
    int start_position;
    int end_position;
    int n_motif;
    std::string sequence;

    // scores per strand
    double ll_unmethylated[2];
    double ll_methylated[2];
    int strands_scored;

    //
    static bool sort_by_position(const ScoredSite& a, const ScoredSite& b) { return a.start_position < b.start_position; }
};

struct MethylationCallingParameters
{
    int min_separation = 10;
    int min_flank = 10;
    std::string methylation_type = "cpg";
    const Alphabet* alphabet = NULL;
};


// Process methylation sites in the provided read
void calculate_methylation_for_read(const OutputHandles& handles,
                                    SquiggleRead& sr,
                                    const MethylationCallingParameters& calling_parameters,
                                    const faidx_t* fai,
                                    const bam_hdr_t* hdr,
                                    const bam1_t* record,
                                    size_t read_idx,
                                    int region_start,
                                    int region_end);

#endif
