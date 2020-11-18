//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_fast5_io -- lightweight functions
// to read specific data from fast5 files
//
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sstream>
#include "nanopolish_fast5_io.h"

//#define DEBUG_FAST5_IO 1

#define LEGACY_FAST5_RAW_ROOT "/Raw/Reads/"

#define H5Z_FILTER_VBZ 32020 //We need to find out what the numerical value for this is

int verbose = 0;

//
fast5_file fast5_open(const std::string& filename)
{
    fast5_file fh;
    fh.hdf5_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Check for attribute that indicates whether it is single or multi-fast5
    // see: https://community.nanoporetech.com/posts/multi-fast5-format
    const std::string indicator_p1 = "/UniqueGlobalKey/";
    const std::string indicator_p2 = indicator_p1 + "tracking_id/";
    bool has_indicator = H5Lexists(fh.hdf5_file, indicator_p1.c_str(), H5P_DEFAULT) && H5Lexists(fh.hdf5_file, indicator_p2.c_str(), H5P_DEFAULT);
    fh.is_multi_fast5 = !has_indicator;
    return fh;
}

//
bool fast5_is_open(fast5_file& fh)
{
    return fh.hdf5_file >= 0;
}

//
void fast5_close(fast5_file& fh)
{
    H5Fclose(fh.hdf5_file);
}

//
std::vector<std::string> fast5_get_multi_read_groups(fast5_file& fh)
{
    std::vector<std::string> out;
    size_t buffer_size = 0;
    char* buffer = NULL;

    // get the number of groups in the root group
    H5G_info_t group_info;
    int ret = H5Gget_info_by_name(fh.hdf5_file, "/", &group_info, H5P_DEFAULT);
    if(ret < 0) {
        fprintf(stderr, "error getting group info\n");
        exit(EXIT_FAILURE);
    }

    for(size_t group_idx = 0; group_idx < group_info.nlinks; ++group_idx) {

        // retrieve the size of this group name
        ssize_t size = H5Lget_name_by_idx(fh.hdf5_file, "/", H5_INDEX_NAME, H5_ITER_INC, group_idx, NULL, 0, H5P_DEFAULT);

        if(size < 0) {
            fprintf(stderr, "error getting group name size\n");
            exit(EXIT_FAILURE);
        }
        size += 1; // for null terminator
           
        if(size > buffer_size) {
            buffer = (char*)realloc(buffer, size);
            buffer_size = size;
        }
    
        // copy the group name
        H5Lget_name_by_idx(fh.hdf5_file, "/", H5_INDEX_NAME, H5_ITER_INC, group_idx, buffer, buffer_size, H5P_DEFAULT);
        buffer[size - 1] = '\0';
        out.push_back(buffer);
    }

    free(buffer);
    buffer = NULL;
    buffer_size = 0;
    return out;
}

//
std::string fast5_get_read_id_single_fast5(fast5_file& fh)
{
    // this function is not supported for multi-fast5 files
    assert(!fh.is_multi_fast5);

    int ret;
    hid_t read_name_attribute, raw_group, attribute_type;
    size_t storage_size = 0;
    char* read_name_str = NULL;

    std::string out = "";
    
    // Get the path to the raw read group
    std::string raw_read_group = fast5_get_raw_read_group(fh, "");
    if(raw_read_group == "") {
        return out;
    }

    return fast5_get_string_attribute(fh, raw_read_group, "read_id");
}

//
raw_table fast5_get_raw_samples(fast5_file& fh, const std::string& read_id, fast5_raw_scaling scaling)
{
    float* rawptr = NULL;
    hid_t space;
    hsize_t nsample;
    herr_t status;
    float raw_unit;
    raw_table rawtbl = { 0, 0, 0, NULL };

    // mostly from scrappie
    std::string raw_read_group = fast5_get_raw_read_group(fh, read_id);

    // Create data set name
    std::string signal_path = raw_read_group + "/Signal";

    hid_t dset = H5Dopen(fh.hdf5_file, signal_path.c_str(), H5P_DEFAULT);
    if (dset < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open dataset '%s' to read raw signal from.\n", signal_path.c_str());
#endif
        goto cleanup2;
    }

    space = H5Dget_space(dset);
    if (space < 0) {
        fprintf(stderr, "Failed to create copy of dataspace for raw signal %s.\n", signal_path.c_str());
        goto cleanup3;
    }

    H5Sget_simple_extent_dims(space, &nsample, NULL);
    rawptr = (float*)calloc(nsample, sizeof(float));
    status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawptr);

    if (status < 0) {
	if(fast5_is_vbz_compressed(fh, read_id) == 1) {
	    fprintf(stderr, "The fast5 file is compressed with VBZ but the required plugin is not loaded. Please read the instructions here: https://github.com/nanoporetech/vbz_compression/issues/5\n");
	    exit(EXIT_FAILURE);
	}
        free(rawptr);
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to read raw data from dataset %s.\n", signal_path.c_str());
#endif
        goto cleanup4;
    }

    // convert to pA
    rawtbl = (raw_table) { nsample, 0, nsample, rawptr };
    raw_unit = scaling.range / scaling.digitisation;
    for (size_t i = 0; i < nsample; i++) {
        rawptr[i] = (rawptr[i] + scaling.offset) * raw_unit;
    }

 cleanup4:
    H5Sclose(space);
 cleanup3:
    H5Dclose(dset);
 cleanup2:
    return rawtbl;
}

std::string fast5_get_context_tags_group(fast5_file& fh, const std::string& read_id)
{
    std::string group = fh.is_multi_fast5 ? "/read_" + read_id + "/context_tags"
                                          : "/UniqueGlobalKey/context_tags";
    return group;
}

//
std::string fast5_get_sequencing_kit(fast5_file& fh, const std::string& read_id)
{
    std::string group = fast5_get_context_tags_group(fh, read_id);
    return fast5_get_string_attribute(fh, group.c_str(), "sequencing_kit");
}

std::string fast5_get_experiment_type(fast5_file& fh, const std::string& read_id)
{
    std::string group = fast5_get_context_tags_group(fh, read_id);
    return fast5_get_string_attribute(fh, group.c_str(), "experiment_type");
}


// from scrappie
float fast5_read_float_attribute(hid_t group, const char *attribute) {
    float val = NAN;
    if (group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Invalid group passed to %s:%d.", __FILE__, __LINE__);
#endif
        return val;
    }

    hid_t attr = H5Aopen(group, attribute, H5P_DEFAULT);
    if (attr < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open attribute '%s' for reading.", attribute);
#endif
        return val;
    }

    H5Aread(attr, H5T_NATIVE_FLOAT, &val);
    H5Aclose(attr);

    return val;
}

uint64_t fast5_read_uint64_attribute(hid_t group, const char *attribute) {
    uint64_t val = 0;
    if (group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Invalid group passed to %s:%d.", __FILE__, __LINE__);
#endif
        return val;
    }

    hid_t attr = H5Aopen(group, attribute, H5P_DEFAULT);
    if (attr < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open attribute '%s' for reading.", attribute);
#endif
        return val;
    }

    H5Aread(attr, H5T_NATIVE_ULLONG, &val);
    H5Aclose(attr);
    return val;
}

uint64_t fast5_get_start_time(fast5_file& fh, const std::string& read_id)
{
    std::string raw_read_group = fast5_get_raw_read_group(fh, read_id);
    hid_t group = H5Gopen(fh.hdf5_file, raw_read_group.c_str(), H5P_DEFAULT);
    if (group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open group %s\n", raw_read_group.c_str());
#endif
        return 0;
    }
    uint64_t t = fast5_read_uint64_attribute(group, "start_time");
    H5Gclose(group);
    return t;
}

//
fast5_raw_scaling fast5_get_channel_params(fast5_file& fh, const std::string& read_id)
{
    // from scrappie
    fast5_raw_scaling scaling = { 0, NAN, NAN, NAN, NAN };

    std::string scaling_path = fh.is_multi_fast5 ? "/read_" + read_id + "/channel_id"
                                                 :  "/UniqueGlobalKey/channel_id";

    hid_t scaling_group = H5Gopen(fh.hdf5_file, scaling_path.c_str(), H5P_DEFAULT);
    if (scaling_group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open group %s\n", scaling_path.c_str());
#endif
        return scaling;
    }
    
    //TODO: open group once?

    // channel
    std::string tmp_id = fast5_get_string_attribute(fh, scaling_path, "channel_number");
    std::stringstream parser(tmp_id);
    parser >> scaling.channel_id;

    scaling.digitisation = fast5_read_float_attribute(scaling_group, "digitisation");
    scaling.offset = fast5_read_float_attribute(scaling_group, "offset");
    scaling.range = fast5_read_float_attribute(scaling_group, "range");
    scaling.sample_rate = fast5_read_float_attribute(scaling_group, "sampling_rate");

    H5Gclose(scaling_group);

    return scaling;
}

//
// Internal functions
//

//
std::string fast5_get_raw_root(fast5_file& fh, const std::string& read_id)
{
    return fh.is_multi_fast5 ? "/read_" + read_id + "/Raw" : "/Raw/Reads";
}

//
std::string fast5_get_raw_read_group(fast5_file& fh, const std::string& read_id)
{
    if(fh.is_multi_fast5) {
        return "/read_" + read_id + "/Raw";
    } else {
        std::string internal_read_name = fast5_get_raw_read_internal_name(fh);
        return internal_read_name != "" ? std::string(LEGACY_FAST5_RAW_ROOT) + "/" + internal_read_name : "";
    }
}

//
std::string fast5_get_raw_read_internal_name(fast5_file& fh)
{
    // This code is From scrappie's fast5_interface

    // retrieve the size of the read name
    ssize_t size =
        H5Lget_name_by_idx(fh.hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, NULL, 0, H5P_DEFAULT);

    if (size < 0) {
        return "";
    }

    // copy the read name out of the fast5
    char* name = (char*)calloc(1 + size, sizeof(char));
    H5Lget_name_by_idx(fh.hdf5_file, LEGACY_FAST5_RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, name, 1 + size, H5P_DEFAULT);

    // cleanup
    std::string out(name);
    free(name);
    return out;
}

//
std::string fast5_get_string_attribute(fast5_file& fh, const std::string& group_name, const std::string& attribute_name)
{
    hid_t group, attribute, attribute_type, native_type;
    std::string out;

    // according to http://hdf-forum.184993.n3.nabble.com/check-if-dataset-exists-td194725.html
    // we should use H5Lexists to check for the existence of a group/dataset using an arbitrary path
    // HDF5 1.8 returns 0 on the root group, so we explicitly check for it
    int ret = group_name == "/" ? 1 : H5Lexists(fh.hdf5_file, group_name.c_str(), H5P_DEFAULT);
    if(ret <= 0) {
        return "";
    }

    // Open the group containing the attribute we want
    group = H5Gopen(fh.hdf5_file, group_name.c_str(), H5P_DEFAULT);
    if(group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "could not open group %s\n", group_name.c_str());
#endif
        goto close_group;
    }

    // Ensure attribute exists
    ret = H5Aexists(group, attribute_name.c_str());
    if(ret <= 0) {
        goto close_group;
    }

    // Open the attribute
    attribute = H5Aopen(group, attribute_name.c_str(), H5P_DEFAULT);
    if(attribute < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "could not open attribute: %s\n", attribute_name.c_str());
#endif
        goto close_attr;
    }

    // Get data type and check it is a fixed-length string
    attribute_type = H5Aget_type(attribute);
    if(attribute_type < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "failed to get attribute type %s\n", attribute_name.c_str());
#endif
        goto close_type;
    }

    if(H5Tget_class(attribute_type) != H5T_STRING) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "attribute %s is not a string\n", attribute_name.c_str());
#endif
        goto close_type;
    }

    native_type = H5Tget_native_type(attribute_type, H5T_DIR_ASCEND);
    if(native_type < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "failed to get native type for %s\n", attribute_name.c_str());
#endif
        goto close_native_type;
    }

    if(H5Tis_variable_str(attribute_type) > 0) {
        // variable length string
        char* buffer;
        ret = H5Aread(attribute, native_type, &buffer);
        if(ret < 0) {
            fprintf(stderr, "error reading attribute %s\n", attribute_name.c_str());
            exit(EXIT_FAILURE);
        }
        out = buffer;
        free(buffer);
        buffer = NULL;

    } else {
        // fixed length string
        size_t storage_size;
        char* buffer;

        // Get the storage size and allocate
        storage_size = H5Aget_storage_size(attribute);
        buffer = (char*)calloc(storage_size + 1, sizeof(char));

        // finally read the attribute
        ret = H5Aread(attribute, attribute_type, buffer);
        if(ret >= 0) {
            out = buffer;
        }

        // clean up
        free(buffer);
    }

close_native_type:
    H5Tclose(native_type);    
close_type:
    H5Tclose(attribute_type);
close_attr:
    H5Aclose(attribute);
close_group:
    H5Gclose(group);

    return out;
}

uint8_t fast5_is_vbz_compressed(fast5_file& fh, const std::string& read_id) {

    hid_t dset, dcpl; 
    H5Z_filter_t filter_id = 0;
    char filter_name[80];
    size_t nelmts = 1; /* number of elements in cd_values */
    unsigned int values_out[1] = {99}; 
    unsigned int flags;

    // mostly from scrappie
    std::string raw_read_group = fast5_get_raw_read_group(fh, read_id);

    // Create data set name
    std::string signal_path = raw_read_group + "/Signal";

    dset = H5Dopen (fh.hdf5_file, signal_path.c_str(), H5P_DEFAULT);

    dcpl = H5Dget_create_plist (dset);

    filter_id = H5Pget_filter2 (dcpl, (unsigned) 0, &flags, &nelmts, values_out, sizeof(filter_name) - 1, filter_name, NULL);

    H5Pclose (dcpl);
    H5Dclose (dset);

    if(filter_id == H5Z_FILTER_VBZ)
        return 1;
    else 
        return 0;
}
