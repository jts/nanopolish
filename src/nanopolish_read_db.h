//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_read_db -- database of reads and their
// associated signal data
//
#ifndef NANOPOLISH_READ_DB
#define NANOPOLISH_READ_DB

#include <map>

struct ReadDBData
{
    // path to the signal-level data for this read
    std::string signal_data_path;
};

class ReadDB
{
    public:
        ReadDB(const std::string& input_reads_filename);
        
        void save() const;
        void load(const std::string& reads_filename);

        void add_raw_signal_path(const std::string& read_id, const std::string& path);
        
        // returns true if all reads in the database have paths to their signal-level data
        bool check_signal_paths() const;
        
        // returns the number of reads in the database
        size_t get_num_reads() const { return m_data.size(); }

        // print some summary stats about the database
        void print_stats() const;

    private:
        
        //
        void import_fastx(const std::string& input_filename, const std::string& output_fasta_filename);

        //
        std::string m_reads_filename;

        //
        std::map<std::string, ReadDBData> m_data;
};

#endif
