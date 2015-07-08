//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_fast5_map - a simple map from a read
// name to a fast5 file path
#ifndef NANOPOLISH_FAST5_MAP
#define NANOPOLISH_FAST5_MAP

#include <string>
#include <map>

class Fast5Map
{
    public:
        Fast5Map(const std::string& fasta_filename);
        
        // return the path for the given read name
        // if the read does not exist in the map, emits an error
        // and exits
        std::string get_path(const std::string& read_name) const;

    private:

        // Read the read -> path map from the header of a fasta file
        void load_from_fasta(std::string fasta_filename);

        // Read the map from a pre-computed .fofn file
        void load_from_fofn(std::string fofn_filename);
        
        // Write the map to the .fofn file
        void write_to_fofn(std::string fofn_filename);

        std::map<std::string, std::string> read_to_path_map;

};

#endif
