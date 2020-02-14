//---------------------------------------------------------
// Copyright 2019 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_fast5_loader -- A class that manages
// opening and reading from fast5 files, in parallel
//
#ifndef NANOPOLISH_FAST5_LOADER_H
#define NANOPOLISH_FAST5_LOADER_H

#include <map>
#include "nanopolish_fast5_loader.h"
#include "nanopolish_fast5_io.h"

// structure holding all of the information we want to get out of a fast5 file
struct Fast5Data
{
    bool is_valid;

    std::string read_name;
    std::string sequencing_kit;
    std::string experiment_type;
    fast5_raw_scaling channel_params;

    // this is allocated and must be freed by the calling function
    raw_table rt;

    uint64_t start_time;
};

class Fast5Loader
{
    public:

        static Fast5Data load_read(const std::string& filename, const std::string& read_name);

        // destructor
        ~Fast5Loader();

    private:

        // singleton accessor function
        static Fast5Loader& getInstance()
        {
            static Fast5Loader instance;
            return instance;
        }

        Fast5Loader();

        // do not allow copies of this classs
        Fast5Loader(Fast5Loader const&) = delete;
        void operator=(Fast5Loader const&) = delete;
};

#endif
