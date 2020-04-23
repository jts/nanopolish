//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_fast5_loader -- A class that manages
// opening and reading from fast5 files
//
#include <omp.h>
#include "nanopolish_fast5_loader.h"

//
Fast5Loader::Fast5Loader()
{

}

//
Fast5Loader::~Fast5Loader()
{

}

Fast5Data Fast5Loader::load_read(const std::string& filename, const std::string& read_name)
{
    Fast5Data data;
    data.rt.n = 0;
    data.rt.raw = NULL;

    fast5_file f5_file = fast5_open(filename);
    if(!fast5_is_open(f5_file)) {
        data.is_valid = false;
        return data;
    }

    data.read_name = read_name;
    data.is_valid = true;

    // metadata
    data.sequencing_kit = fast5_get_sequencing_kit(f5_file, read_name);
    data.experiment_type = fast5_get_experiment_type(f5_file, read_name);

    // raw data
    data.channel_params = fast5_get_channel_params(f5_file, read_name);
    data.rt = fast5_get_raw_samples(f5_file, read_name, data.channel_params);

    data.start_time = fast5_get_start_time(f5_file, read_name);
    fast5_close(f5_file);

    return data;
}
