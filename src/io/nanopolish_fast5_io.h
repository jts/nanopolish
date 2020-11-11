//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_fast5_io -- lightweight functions
// to read specific data from fast5 files. Some functions
// ported from scrappie's fast5_interface.c
//
#ifndef NANOPOLISH_FAST5_IO_H
#define NANOPOLISH_FAST5_IO_H

#include <string>
#include <vector>
#include <hdf5.h>

extern "C" {
#include "event_detection.h"
#include "scrappie_common.h"
}

// From scrappie
typedef struct {

    int channel_id;
    // Parameters for scaling raw data from ADC values to pA
    float digitisation;
    float offset;
    float range;
    float sample_rate;
} fast5_raw_scaling;

//
struct fast5_file
{
    hid_t hdf5_file;
    bool is_multi_fast5;
};

//
// External API
//

//
// Basic I/O
//

// open the file and return handle
fast5_file fast5_open(const std::string& filename);

// check if file is open
bool fast5_is_open(fast5_file& fh);

// close the file
void fast5_close(fast5_file& fh);

//
// Functions to get the names of the read(s) contained in the file
//

// Get the identifier of a read from the hdf5 file (eg 00041f-....)
std::string fast5_get_read_id_single_fast5(fast5_file& fh);

// Get a vector of read groups for a multi-fast5 file (eg [read_00041f-..., read_1243fe-....])
std::vector<std::string> fast5_get_multi_read_groups(fast5_file& fh);

//
// Functions to get the samples or metadata
//

// get the raw samples from this file
raw_table fast5_get_raw_samples(fast5_file& fh, const std::string& read_id, fast5_raw_scaling scaling);

// Get the sequencing kit
std::string fast5_get_sequencing_kit(fast5_file& fh, const std::string& read_id);

// Get the experiment type attribute
std::string fast5_get_experiment_type(fast5_file& fh, const std::string& read_id);

// Get sample rate, and ADC-to-pA scalings
fast5_raw_scaling fast5_get_channel_params(fast5_file& fh, const std::string& read_id);

// Get the start time of this read
uint64_t fast5_get_start_time(fast5_file& fh, const std::string& read_id);

//
// Internal utility functions
//

// get the name of the raw read in the file (eg Read_1234)
std::string fast5_get_raw_read_internal_name(fast5_file& fh);

// get the name of the raw read group (eg /Raw/Read/Read_1234 or /read_00041f-.../Raw/)
std::string fast5_get_raw_read_group(fast5_file& fh, const std::string& read_id);

//
std::string fast5_get_string_attribute(fast5_file& fh, const std::string& group_name, const std::string& attribute_name);

uint8_t fast5_is_vbz_compressed(fast5_file& fh, const std::string& read_id);

#endif
