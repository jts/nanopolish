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
    //  Information for scaling raw data from ADC values to pA
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

// open the file and return the hdf ID
fast5_file fast5_open(const std::string& filename);

bool fast5_is_open(fast5_file& fh);

// close the file
void fast5_close(fast5_file& fh);

// get the raw samples from this file
raw_table fast5_get_raw_samples(fast5_file& fh, fast5_raw_scaling scaling);

// get the name of the raw read in the file (eg Read_1234)
std::string fast5_get_raw_read_name(fast5_file& fh);

// get the name of the raw read group (eg /Raw/Read/Read_1234)
std::string fast5_get_raw_read_group(fast5_file& fh);

// Get the identifier of a read from the hdf5 file
std::string fast5_get_read_id(fast5_file& fh);

// Get the experiment type attribute
std::string fast5_get_experiment_type(fast5_file& fh);

// Get sample rate, and ADC-to-pA scalings
fast5_raw_scaling fast5_get_channel_params(fast5_file& fh);

//
// Internal utility functions
//
std::string fast5_get_string_attribute(fast5_file& fh, const std::string& group_name, const std::string& attribute_name);

#endif
