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

Fast5Data Fast5Loader::load_read(slow5_file_t *slow5_file, const std::string &read_name)
{
    Fast5Data data;
    slow5_rec_t *rec = NULL;
    int ret = slow5_get(read_name.c_str(), &rec, slow5_file);
    if(ret < 0){
        fprintf(stderr,"Error in when fetching the read\n");
        exit(EXIT_FAILURE);
    }

    data.rt.n = rec->len_raw_signal;
    data.channel_params.range = rec->range;
    data.channel_params.offset = rec->offset;
    data.channel_params.digitisation = rec->digitisation;
    data.channel_params.sample_rate = rec->sampling_rate;

    int err;
    char *cid = slow5_aux_get_string(rec, "channel_number", NULL, &err);
    if(err < 0){
        fprintf(stderr,"[warning] Error in when fetching the channel_number\n");
    }else{
        data.channel_params.channel_id = atoi(cid);
    }
    data.start_time = slow5_aux_get_uint64(rec, "start_time", &err);
    if(err < 0){
        fprintf(stderr,"[warning] Error in when fetching the start_time\n");
    }
    data.read_name = read_name;
    // metadata
    char* sequencing_kit = slow5_hdr_get("sequencing_kit", 0, slow5_file->header);
    data.sequencing_kit = (sequencing_kit)?std::string(sequencing_kit):"";
    char* experiment_type = slow5_hdr_get("experiment_type", 0, slow5_file->header);
    data.experiment_type = (experiment_type)?std::string(experiment_type):"";

    // raw data
    // convert to pA
    float* rawptr = (float*)calloc(rec->len_raw_signal, sizeof(float));
    raw_table rawtbl = { 0, 0, 0, NULL };
    data.rt = (raw_table) { rec->len_raw_signal, 0, rec->len_raw_signal, rawptr };
    hsize_t nsample = rec->len_raw_signal;
    float digitisation = rec->digitisation;
    float offset = rec->offset;
    float range = rec->range;
    float raw_unit = range / digitisation;
    for (size_t i = 0; i < nsample; i++) {
        float signal = rec->raw_signal[i];
        rawptr[i] = (signal + offset) * raw_unit;
    }
    slow5_rec_free(rec);
    data.is_valid = true;
    return data;
}
