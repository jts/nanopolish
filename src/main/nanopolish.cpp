//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish.cpp -- main driver program
//
#include <string>
#include <map>
#include <functional>
#include "logsum.h"
#include "nanopolish_index.h"
#include "nanopolish_extract.h"
#include "nanopolish_call_variants.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_getmodel.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_call_methylation.h"
#include "nanopolish_scorereads.h"
#include "nanopolish_phase_reads.h"
#include "nanopolish_train_poremodel_from_basecalls.h"

int print_usage(int argc, char **argv);
int print_version(int argc, char **argv);

static std::map< std::string, std::function<int(int, char**)> > programs = {
    {"help",        print_usage},
    {"--help",      print_usage},
    {"--version",   print_version},
    {"index",       index_main},
    {"extract",     extract_main},
    {"eventalign",  eventalign_main},
    {"getmodel",    getmodel_main},
    {"variants",    call_variants_main},
    {"methyltrain", methyltrain_main},
    {"scorereads",  scorereads_main} ,
    {"phase-reads",  phase_reads_main} ,
    {"call-methylation",  call_methylation_main}
};

int print_usage(int, char **)
{
    std::cout << "usage: nanopolish [command] [options]" << std::endl;
    std::cout << "  valid commands: " << std::endl;
    for (const auto &item : programs){
        std::cout << "    " << item.first << std::endl;
    }
    std::cout << "  for help on given command, type nanopolish command --help" << std::endl;
    return 0;
}

int print_version(int, char **)
{
    static const char *VERSION_MESSAGE =
    "nanopolish version " PACKAGE_VERSION "\n"
    "Written by Jared Simpson.\n"
    "\n"
    "Copyright 2015-2017 Ontario Institute for Cancer Research\n";
    std::cout << VERSION_MESSAGE << std::endl;
    return 0;
}

int main(int argc, char** argv)
{
    // Turn off HDF's exception printing, which is generally unhelpful for users
    H5Eset_auto(0, NULL, NULL);

    int ret = 0;
    if(argc <= 1) {
        printf("error: no command provided\n");
        print_usage(argc - 1 , argv + 1);
        return 0;
    } else {
        std::string command(argv[1]);
        auto iter = programs.find(command);
        if (iter != programs.end()) 
            ret = iter->second( argc - 1, argv + 1);
        else
            ret = print_usage( argc - 1, argv + 1);
    }


    // Emit a warning when some reads had to be skipped
    extern int g_total_reads;
    extern int g_unparseable_reads;
    extern int g_qc_fail_reads;
    extern int g_failed_calibration_reads;
    extern int g_failed_alignment_reads;
    extern int g_bad_fast5_file;
    if(g_total_reads > 0) {
        fprintf(stderr, "[post-run summary] total reads: %d, unparseable: %d, qc fail: %d, could not calibrate: %d, no alignment: %d, bad fast5: %d\n", 
            g_total_reads, g_unparseable_reads, g_qc_fail_reads, g_failed_calibration_reads, g_failed_alignment_reads, g_bad_fast5_file);
    }
    return ret;
}
