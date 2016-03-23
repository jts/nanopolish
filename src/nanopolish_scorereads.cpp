//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_scorereads -- score reads against an alignment, model
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include <cstddef>
#include "htslib/faidx.h"
#include "nanopolish_alphabet.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_iupac.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_anchor.h"
#include "nanopolish_fast5_map.h"
#include "H5pubconf.h"

//
// Getopt
//
#define SUBPROGRAM "scorereads"

static const char *SCOREREADS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *SCOREREADS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] --reads reads.fa --bam alignments.bam --genome genome.fa\n"
"Score reads against an alignment, model\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --version                        display version\n"
"      --help                           display this help and exit\n"
"  -m, --models-fofn=FILE               optionally use these models rather than models in fast5\n"
"  -c  --calibrate                      recalibrate aligned reads to model before scoring\n"
"  -i  --individual-reads=READ,READ     optional comma-delimited list of readnames to score\n"
"  -r, --reads=FILE                     the 2D ONT reads are in fasta FILE\n"
"  -b, --bam=FILE                       the reads aligned to the genome assembly are in bam FILE\n"
"  -g, --genome=FILE                    the genome we are computing a consensus for is in FILE\n"
"  -t, --threads=NUM                    use NUM threads (default: 1)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static unsigned int calibrate=0;
    static std::string reads_file;
    static std::string bam_file;
    static std::string genome_file;
    static std::string models_fofn;
    static std::string region;
    static std::vector<std::string> readnames;
    static int num_threads = 1;
    static int batch_size = 128;
}

static const char* shortopts = "i:r:b:g:t:m:vc";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",            no_argument,       NULL, 'v' },
    { "calibrate",          no_argument,       NULL, 'c' },
    { "reads",              required_argument, NULL, 'r' },
    { "bam",                required_argument, NULL, 'b' },
    { "genome",             required_argument, NULL, 'g' },
    { "threads",            required_argument, NULL, 't' },
    { "models-fofn",        required_argument, NULL, 'm' },
    { "individual-reads",   required_argument, NULL, 'i' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

double model_score(SquiggleRead &sr,
                   const size_t strand_idx,
                   const faidx_t *fai, 
                   const std::vector<EventAlignment> &alignment_output,
                   const size_t events_per_segment)  
{
    double curr_score = 0;
    size_t nevents = 0;

    for(int align_start_idx = events_per_segment; 
               align_start_idx < (int)alignment_output.size() - (int)events_per_segment; 
               align_start_idx += events_per_segment) {

        const EventAlignment& align_start = alignment_output[align_start_idx];
        const EventAlignment& align_end = alignment_output[align_start_idx + events_per_segment];
        std::string contig = alignment_output.front().ref_name.c_str();

        // Set up event data
        HMMInputData data;
        data.read = &sr;
        data.anchor_index = -1; // unused
        data.strand = strand_idx;
        data.rc = alignment_output.front().rc;
        data.event_start_idx = align_start.event_idx;
        data.event_stop_idx = align_end.event_idx;
        data.event_stride = data.event_start_idx <= data.event_stop_idx ? 1 : -1;
        
        // Set up reference data
        int ref_start_pos = align_start.ref_position;
        int ref_end_pos = align_end.ref_position;
        int fetched_len = 0;

        assert(ref_end_pos >= ref_start_pos);

        // Extract the reference sequence for this region
        std::string ref_seq = get_reference_region_ts(fai, contig.c_str(), ref_start_pos, 
                                                      ref_end_pos, &fetched_len);

        if (fetched_len <= (int)sr.pore_model[strand_idx].k)
            continue;

        const Alphabet *alphabet = sr.pore_model[strand_idx].pmalphabet;
    
        ref_seq = alphabet->disambiguate(ref_seq);
        HMMInputSequence sequence(ref_seq, alphabet->reverse_complement(ref_seq), alphabet);

        // Run HMM using current model
        curr_score += profile_hmm_score(sequence, data, 0);
        nevents += events_per_segment;
    }

    if (nevents == 0)
        return +1;
    else
        return curr_score/nevents;
}

std::vector<EventAlignment> alignment_from_read(SquiggleRead& sr,
                                                const size_t strand_idx, 
                                                const size_t read_idx,
                                                const ModelMap& model_map,
                                                const faidx_t* fai,                        
                                                const bam_hdr_t* hdr,                        
                                                const bam1_t* record, 
                                                int region_start,
                                                int region_end)
{
    // optionally replace model  
    if (!opt::models_fofn.empty()) {
        std::string curr_model = sr.pore_model[strand_idx].name;
        auto model_iter = model_map.find(curr_model);

        if(model_iter != model_map.end()) {
            sr.pore_model[strand_idx].update_states(model_iter->second);
        } else {
            printf("Error: model %s not found\n", curr_model.c_str());
            assert(false && "Model not found");
        }
    }        

    // set k
    //uint32_t k = sr.pore_model[strand_idx].k;

    // Align to the new model
    EventAlignmentParameters params;
    params.sr = &sr;
    params.fai = fai;
    params.hdr = hdr;
    params.record = record;
    params.strand_idx = strand_idx;

    params.alphabet = sr.pore_model[strand_idx].pmalphabet;
    params.read_idx = read_idx;
    params.region_start = region_start;
    params.region_end = region_end;
    return align_read_to_ref(params);
} 

void parse_scorereads_options(int argc, char** argv)
{
    bool die = false;
    std::string readlist;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case 'r': arg >> opt::reads_file; break;
            case 'b': arg >> opt::bam_file; break;
            case 'g': arg >> opt::genome_file; break;
            case 't': arg >> opt::num_threads; break;
            case 'm': arg >> opt::models_fofn; break;
            case 'i': arg >> readlist; break;
            case 'v': opt::verbose++; break;
            case 'c': opt::calibrate = 1; break;
            case '?': die = true; break;
            case OPT_HELP:
                std::cout << SCOREREADS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << SCOREREADS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind > 0) {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if(opt::num_threads <= 0) {
        std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::num_threads << "\n";
        die = true;
    }

    if(opt::reads_file.empty()) {
        std::cerr << SUBPROGRAM ": a --reads file must be provided\n";
        die = true;
    }
    
    if(opt::genome_file.empty()) {
        std::cerr << SUBPROGRAM ": a --genome file must be provided\n";
        die = true;
    }

    if(opt::bam_file.empty()) {
        std::cerr << SUBPROGRAM ": a --bam file must be provided\n";
        die = true;
    }
    
    // this is much cleaner with sregex_token_iterator, which isn't implemented in gcc until 4.9
    if (!readlist.empty()) {
        size_t start = readlist.find_first_not_of(","), end=start;
        while (start != std::string::npos){
                end = readlist.find(",", start);
                opt::readnames.push_back(readlist.substr(start, end-start));
                start = readlist.find_first_not_of(",", end);
        }
    }

    if (die) 
    {
        std::cout << "\n" << SCOREREADS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}


int scorereads_main(int argc, char** argv)
{
    parse_scorereads_options(argc, argv);
    omp_set_num_threads(opt::num_threads);

    Fast5Map name_map(opt::reads_file);
    ModelMap models;
    if (!opt::models_fofn.empty())
        models = read_models_fofn(opt::models_fofn);
    
    // Open the BAM and iterate over reads

    // load bam file
    htsFile* bam_fh = sam_open(opt::bam_file.c_str(), "r");
    assert(bam_fh != NULL);

    // load bam index file
    std::string index_filename = opt::bam_file + ".bai";
    hts_idx_t* bam_idx = bam_index_load(index_filename.c_str());
    assert(bam_idx != NULL);

    // read the bam header
    bam_hdr_t* hdr = sam_hdr_read(bam_fh);
    
    // load reference fai file
    faidx_t *fai = fai_load(opt::genome_file.c_str());

    hts_itr_t* itr;

    // If processing a region of the genome, only emit events aligned to this window
    int clip_start = -1;
    int clip_end = -1;

    if(opt::region.empty()) {
        // TODO: is this valid?
        itr = sam_itr_queryi(bam_idx, HTS_IDX_START, 0, 0);
    } else {

        fprintf(stderr, "Region: %s\n", opt::region.c_str());
        itr = sam_itr_querys(bam_idx, hdr, opt::region.c_str());
        hts_parse_reg(opt::region.c_str(), &clip_start, &clip_end);
    }

#ifndef H5_HAVE_THREADSAFE
    if(opt::num_threads > 1) {
        fprintf(stderr, "You enabled multi-threading but you do not have a threadsafe HDF5\n");
        fprintf(stderr, "Please recompile nanopolish's built-in libhdf5 or run with -t 1\n");
        exit(1);
    }
#endif

    // Initialize iteration
    std::vector<bam1_t*> records(opt::batch_size, NULL);
    for(size_t i = 0; i < records.size(); ++i) {
        records[i] = bam_init1();
    }

    int result;
    size_t num_reads_realigned = 0;
    size_t num_records_buffered = 0;

    do {
        assert(num_records_buffered < records.size());
        
        // read a record into the next slot in the buffer
        result = sam_itr_next(bam_fh, itr, records[num_records_buffered]);
        num_records_buffered += result >= 0;

        // realign if we've hit the max buffer size or reached the end of file
        if(num_records_buffered == records.size() || result < 0) {
            #pragma omp parallel for schedule(dynamic)
            for(size_t i = 0; i < num_records_buffered; ++i) {
                bam1_t* record = records[i];
                size_t read_idx = num_reads_realigned + i;
                if( (record->core.flag & BAM_FUNMAP) == 0) {

                    //load read
                    std::string read_name = bam_get_qname(record);
                    std::string fast5_path = name_map.get_path(read_name);
                    SquiggleRead sr(read_name, fast5_path);

                    // TODO: early exit when have processed all of the reads in readnames
                    if (!opt::readnames.empty() && 
                         std::find(opt::readnames.begin(), opt::readnames.end(), read_name) == opt::readnames.end() )
                            continue;

                    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {
                        std::vector<EventAlignment> ao = alignment_from_read(sr, strand_idx, read_idx,
                                                                             models, fai, hdr,
                                                                             record, clip_start, clip_end);
                        if (ao.size() == 0)
                            continue;

                        // Update pore model based on alignment
                        if ( opt::calibrate ) 
                            recalibrate_model(sr, strand_idx, ao, false);

                        double score = model_score(sr, strand_idx, fai, ao, 500);
                        if (score > 0) 
                            continue;
                        #pragma omp critical(print)
                        std::cout << read_name << " " << ( strand_idx ? "complement" : "template" ) 
                                  << " " << sr.pore_model[strand_idx].name << " " << score << std::endl;
                    } 
                }
            }

            num_reads_realigned += num_records_buffered;
            num_records_buffered = 0;
        }

    } while(result >= 0);
    
    // cleanup records
    for(size_t i = 0; i < records.size(); ++i) {
        bam_destroy1(records[i]);
    }

    // cleanup
    sam_itr_destroy(itr);
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    sam_close(bam_fh);
    hts_idx_destroy(bam_idx);
    return 0;
}

