//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_poremodel -- Representation of the Oxford
// Nanopore sequencing model, as described in a FAST5 file
//
#include "nanopolish_poremodel.h"
#include <fstream>
#include <sstream>
#include <cstring>
#include <bits/stl_algo.h>
#include "../fast5/src/fast5.hpp"

void PoreModel::bake_gaussian_parameters()
{
    scaled_params.resize(states.size());
    scaled_states.resize(states.size());

    for(int i = 0; i < states.size(); ++i) {

        // calculate the derived sd_lambda parameter
        states[i].sd_lambda = pow(states[i].sd_mean, 3.0) / pow(states[i].sd_stdv, 2.0);

        // as per ONT documents
        scaled_states[i].level_mean = states[i].level_mean * scale + shift;
        scaled_states[i].level_stdv = states[i].level_stdv * var;

        scaled_states[i].sd_mean = states[i].sd_mean * scale_sd;
        scaled_states[i].sd_lambda = states[i].sd_lambda * var_sd;
        scaled_states[i].sd_stdv = sqrt(pow(scaled_states[i].sd_mean, 3) / scaled_states[i].sd_lambda);

        // for efficiency
        scaled_states[i].level_log_stdv = log(scaled_states[i].level_stdv);
        scaled_states[i].sd_log_lambda = log(scaled_states[i].sd_lambda);

        // for compatibility
        scaled_params[i].mean = scaled_states[i].level_mean;
        scaled_params[i].stdv = scaled_states[i].level_stdv;
        scaled_params[i].log_stdv = scaled_states[i].level_log_stdv;
    }
    is_scaled = true;
}

void add_found_bases(char *known, const char *kmer) {
    char newbase[2];
    int posn;
    newbase[1] = '\0';

    while ( (posn = strspn(kmer, known)) != strlen(kmer) ){
        newbase[0] = kmer[posn];
        strcat(known, newbase);
    }
    return;
}

PoreModel::PoreModel(const std::string filename, const Alphabet *alphabet) : pmalphabet(alphabet), is_scaled(false)
{
    model_filename = filename;
    std::ifstream model_reader(filename);
    std::string model_line;

    bool firstKmer = true;
    int ninserted = 0;

    shift_offset = 0.0f;

    const size_t maxNucleotides=50;
    char bases[maxNucleotides+1]="";

    std::map<std::string, PoreModelStateParams> kmers;

    while (getline(model_reader, model_line)) {
        std::stringstream parser(model_line);

        // Extract the model name from the header
        if (model_line.find("#model_name") != std::string::npos) {
            std::string dummy;
            parser >> dummy >> name;
        }

        // Extract shift offset from the header
        // This will be applied to the per-read shift values
        // to allow switching between models with different averages
        if (model_line.find("#shift_offset") != std::string::npos) {
            std::string dummy;
            parser >> dummy >> shift_offset;
        }

        // Use the alphabet defined in the header if available
        if (model_line.find("#alphabet") != std::string::npos) {
            std::string dummy;
            std::string alphabet_name;
            parser >> dummy >> alphabet_name;
            pmalphabet = get_alphabet_by_name(alphabet_name);
        }


        // skip the rest of the header
        if (model_line[0] == '#' || model_line.find("kmer") == 0) {
            continue;
        }

        std::string kmer;
        PoreModelStateParams params;
        parser >> kmer >> params.level_mean >> params.level_stdv >> params.sd_mean >> params.sd_stdv;

        kmers[kmer] = params;
        add_found_bases(bases, kmer.c_str());

        if (firstKmer) {
            k = kmer.length();
            firstKmer = false;
        }
    }

    if (pmalphabet == nullptr) 
        pmalphabet = best_alphabet(bases);

    assert( pmalphabet != nullptr );

    states.resize(pmalphabet->get_num_strings(k));
    for (const auto &iter : kmers ) {
        ninserted++;
        states[ pmalphabet->kmer_rank(iter.first.c_str(), k) ] = iter.second;
    }
    assert( ninserted == states.size() );

    is_scaled = false;
}

PoreModel::PoreModel(fast5::File *f_p, const size_t strand, const Alphabet *alphabet) : pmalphabet(alphabet)
{
    const size_t maxNucleotides=50;
    char bases[maxNucleotides+1]="";

    std::map<std::string, PoreModelStateParams> kmers;

    std::vector<fast5::Model_Entry> model = f_p->get_model(strand);
    k = (uint32_t) strlen(model[0].kmer);

    // Copy into the pore model for this read
    for(size_t mi = 0; mi < model.size(); ++mi) {
        const fast5::Model_Entry& curr = model[mi];

        std::string stringkmer(curr.kmer);
        kmers[stringkmer] = { static_cast<float>(curr.level_mean),
                              static_cast<float>(curr.level_stdv),
                              static_cast<float>(curr.sd_mean),
                              static_cast<float>(curr.sd_stdv) };
        add_found_bases(bases, curr.kmer);
    }

    if (pmalphabet == nullptr)
        pmalphabet = best_alphabet(bases);
    assert( pmalphabet != nullptr );

    states.resize( pmalphabet->get_num_strings(k) );
    assert(states.size() == model.size());

    for (const auto &iter : kmers ) {
        states[ pmalphabet->kmer_rank(iter.first.c_str(), k) ] = iter.second;
    }

    // Load the scaling parameters for the pore model
    fast5::Model_Parameters params = f_p->get_model_parameters(strand);
    drift = params.drift;
    scale = params.scale;
    scale_sd = params.scale_sd;
    shift = params.shift;
    var = params.var;
    var_sd = params.var_sd;

    // no offset needed when loading directly from the fast5
    shift_offset = 0.0f;

    // apply shift/scale transformation to the pore model states
    bake_gaussian_parameters();

    // Read and shorten the model name
    std::string temp_name = f_p->get_model_file(strand);
    std::string leader = "/opt/chimaera/model/";

    size_t lp = temp_name.find(leader);
    // leader not found
    if(lp == std::string::npos) {
        name = temp_name;
    } else {
        name = temp_name.substr(leader.size());
    }

    std::replace(name.begin(), name.end(), '/', '_');
}

void PoreModel::write(const std::string filename, const std::string modelname) 
{
    std::string outmodelname = modelname;
    if(modelname.empty())
        outmodelname = name;

    std::ofstream writer(filename);
    writer << "#model_name\t" << outmodelname << std::endl;
    writer << "#shift_offset\t" << shift_offset << std::endl;

    std::string curr_kmer(k,pmalphabet->base(0));
    for(size_t ki = 0; ki < states.size(); ++ki) {
        writer << curr_kmer << "\t" << states[ki].level_mean << "\t" << states[ki].level_stdv << "\t"
               << states[ki].sd_mean << "\t" << states[ki].sd_stdv << std::endl;
        pmalphabet->lexicographic_next(curr_kmer);
    }
    writer.close();
}

void PoreModel::update_states( const PoreModel &other ) 
{
    k = other.k;
    pmalphabet = other.pmalphabet;
    shift += other.shift_offset;
    update_states( other.states );
}

void PoreModel::update_states( const std::vector<PoreModelStateParams> &otherstates ) 
{
    states = otherstates;
    if (is_scaled) {
        bake_gaussian_parameters();
    }
}

ModelMap read_models_fofn(const std::string& fofn_name, const Alphabet *alphabet)
{
    ModelMap out;
    std::ifstream fofn_reader(fofn_name);
	if(!fofn_reader.is_open()) {
		fprintf(stderr, "Error: could not read %s\n", fofn_name.c_str());
		exit(EXIT_FAILURE);
	}
    std::string model_filename;

    while(getline(fofn_reader, model_filename)) {
        printf("reading %s\n", model_filename.c_str());
        PoreModel p(model_filename, alphabet);
        assert(!p.name.empty());
        out[p.name] = p;
    }
    return out;
}

