//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_fast5_io -- lightweight functions
// to read specific data from fast5 files
//
#include <string.h>
#include "nanopolish_fast5_io.h"

#define RAW_ROOT "/Raw/Reads/"

//
hid_t fast5_open(const std::string& filename)
{
    hid_t hdf5file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    return hdf5file;
}

void fast5_close(hid_t hdf5_file)
{
    H5Fclose(hdf5_file);
}

//
std::string fast5_get_raw_read_name(hid_t hdf5_file)
{
    // This code is From scrappie's fast5_interface

    // retrieve the size of the read name
    ssize_t size =
        H5Lget_name_by_idx(hdf5_file, RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, NULL,
                           0, H5P_DEFAULT);

    if (size < 0) {
        return "";
    }

    // copy the read name out of the fast5
    char* name = (char*)calloc(1 + size, sizeof(char));
    H5Lget_name_by_idx(hdf5_file, RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, name,
                       1 + size, H5P_DEFAULT);

    // cleanup
    std::string out(name);
    free(name);
    return out;
}

//
std::string fast5_get_read_id(hid_t hdf5_file)
{
    int ret;
    hid_t read_name_attribute, raw_group, attribute_type;
    size_t storage_size = 0;
    char* read_name_str = NULL;

    std::string out = "";
    
    // Get the path to the raw read
    std::string raw_read = fast5_get_raw_read_name(hdf5_file);
    if(raw_read == "") {
        return out;
    }
    std::string raw_group_path = RAW_ROOT + raw_read;

    // Open the group /Raw/Reads/Read_nnn
    raw_group = H5Gopen(hdf5_file, raw_group_path.c_str(), H5P_DEFAULT);
    if(raw_group < 0) {
        goto close_group;
    }

    // Open the attribute
    read_name_attribute = H5Aopen(raw_group, "read_id", H5P_DEFAULT);
    if(read_name_attribute < 0) {
        goto close_attr;
    }

    // Read the attribute
    storage_size = H5Aget_storage_size(read_name_attribute);
    read_name_str = (char*)calloc(storage_size + 1, sizeof(char));
    attribute_type = H5Aget_type(read_name_attribute);

    ret = H5Aread(read_name_attribute, attribute_type, read_name_str);
    if(ret >= 0) {
        out = read_name_str;
    }

    // clean up
    free(read_name_str);
    H5Tclose(attribute_type);
close_attr:
    H5Aclose(read_name_attribute);
close_group:
    H5Gclose(raw_group);
    return out;
}
