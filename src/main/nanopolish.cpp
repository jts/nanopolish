//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish.cpp -- main driver program
//
#include <stdio.h>
#include <string>
#include "logsum.h"
#include "nanopolish_call_variants.h"
#include "nanopolish_consensus.h"
#include "nanopolish_eventalign.h"
#include "nanopolish_getmodel.h"

// This code needs to be run before any of the program logic
// It sets up pre-computed values and caches
void initialize()
{
    p7_FLogsumInit();
}

void print_usage()
{
    printf("usage: nanopolish [command] [options]\n");
}

int main(int argc, char** argv)
{
    initialize();

    if(argc <= 1) {
        printf("error: no command provided\n");
        print_usage();
        return 0;
    } else {
        std::string command(argv[1]);
        if(command == "help" || command == "--help") {
            print_usage();
            return 0;
        } else if(command == "consensus") {
            consensus_main(argc - 1, argv + 1);
            return 0;
        } else if(command == "eventalign") {
            eventalign_main(argc - 1, argv + 1);
            return 0;
        } else if(command == "getmodel") {
            getmodel_main(argc - 1, argv + 1);
            return 0;
        } else if(command == "variants") {
            call_variants_main(argc - 1, argv + 1);
            return 0;
        }
    }
}
