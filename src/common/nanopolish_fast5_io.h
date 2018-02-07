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

// open the file and return the hdf ID
hid_t fast5_open(const std::string& filename);

// close the file
void fast5_close(hid_t hdf5_file);

// get the name of the raw read in the file (eg Read_1234)
std::string fast5_get_raw_read_name(hid_t hdf5_file);

// Get the identifier of a read from the hdf5 file
std::string fast5_get_read_id(hid_t hdf5_file);

#endif
