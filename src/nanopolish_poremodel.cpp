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
#include <fast5.hpp>

void PoreModel::bake_gaussian_parameters()
{
    scaled_params.resize(states.size());
    scaled_states.resize(states.size());

    for(unsigned i = 0; i < states.size(); ++i) {

        // as per ONT documents
        scaled_states[i].level_mean = states[i].level_mean * scale + shift;
        scaled_states[i].level_stdv = states[i].level_stdv * var;
        scaled_states[i].sd_mean = states[i].sd_mean * scale_sd;
        scaled_states[i].sd_lambda = states[i].sd_lambda * var_sd;
        scaled_states[i].update_sd_stdv();

        // for efficiency
        scaled_states[i].update_logs();

        // for compatibility
        scaled_params[i].mean = scaled_states[i].level_mean;
        scaled_params[i].stdv = scaled_states[i].level_stdv;
        scaled_params[i].log_stdv = scaled_states[i].level_log_stdv;
    }
    is_scaled = true;
}

void add_found_bases(char *known, const char *kmer) {
    char newbase[2];
    unsigned posn;
    newbase[1] = '\0';

    while ( (posn = strspn(kmer, known)) != strlen(kmer) ){
        newbase[0] = kmer[posn];
        strcat(known, newbase);
    }
    return;
}

PoreModel::PoreModel(const std::string filename, const Alphabet *alphabet) : is_scaled(false), pmalphabet(alphabet)
{
    model_filename = filename;
    std::ifstream model_reader(filename);
    std::string model_line;

    bool model_metadata_in_header = false;
    bool firstKmer = true;
    std::string in_kit;
    std::string in_strand;

    unsigned ninserted = 0;

    this->shift = 0.0;
    this->scale = 1.0;
    this->drift = 0.0;
    this->var = 1.0;
    this->scale_sd = 1.0;
    this->var_sd = 1.0;
    this->shift_offset = 0.0f;
    this->scale_offset = 0.0f;

    const size_t maxNucleotides = 50;
    char bases[maxNucleotides+1] = "";

    std::map<std::string, PoreModelStateParams> kmers;
    while (getline(model_reader, model_line)) {
        std::stringstream parser(model_line);

        // Extract the model name from the header
        if (model_line.find("#model_name") != std::string::npos) {
            std::string dummy;
            parser >> dummy >> this->name;
        }

        // Extract the strand from the header
        if (model_line.find("#strand") != std::string::npos) {
            std::string dummy;
            parser >> dummy >> in_strand;
            model_metadata_in_header = true;
        }

        // Extract the sequencing kit version from the header
        if (model_line.find("#kit") != std::string::npos) {
            std::string dummy;
            parser >> dummy >> in_kit;
        }

        if (model_line.find("#type") != std::string::npos) {
            std::string dummy;
            parser >> dummy >> this->type;
        }

        // Extract shift/scale offset from the header
        // This will be applied to the per-read shift values
        // to allow switching between models with different averages
        if (model_line.find("#shift_offset") != std::string::npos) {
            std::string dummy;
            parser >> dummy >> this->shift_offset;
        }
        
        if (model_line.find("#scale_offset") != std::string::npos) {
            std::string dummy;
            parser >> dummy >> this->scale_offset;
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

        // ig_lambda (R9), weight currently not read
        parser >> kmer >> params.level_mean >> params.level_stdv >> params.sd_mean >> params.sd_stdv;

        params.update_sd_lambda();
        params.update_logs();

        kmers[kmer] = params;
        add_found_bases(bases, kmer.c_str());

        if (firstKmer) {
            this->k = kmer.length();
            firstKmer = false;
        }
    }

    if(!model_metadata_in_header) {
        this->metadata = get_model_metadata_from_name(this->name);
    } else {
        set_metadata(in_kit, in_strand);
    }

    assert(metadata.model_idx < 3);

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

PoreModel::PoreModel(fast5::File *f_p, const size_t strand, const std::string& bc_gr, const Alphabet *alphabet) : pmalphabet(alphabet)
{
    const size_t maxNucleotides=50;
    char bases[maxNucleotides+1]="";

    std::map<std::string, PoreModelStateParams> kmers;

    std::vector<fast5::Model_Entry> model = f_p->get_basecall_model(strand, bc_gr);
    k = array2str(model[0].kmer).length();
    assert(k == 5 || k == 6);

    // Copy into the pore model for this read
    for(size_t mi = 0; mi < model.size(); ++mi) {
        const fast5::Model_Entry& curr = model[mi];

        std::string stringkmer = array2str(curr.kmer);
        assert(stringkmer.size() == k);
        kmers[stringkmer] = curr;
        add_found_bases(bases, stringkmer.c_str());
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
    fast5::Model_Parameters params = f_p->get_basecall_model_params(strand, bc_gr);
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
    std::string temp_name = f_p->get_basecall_model_file(strand, bc_gr);
    std::string leader = "/opt/chimaera/model/";

    size_t lp = temp_name.find(leader);
    // leader not found
    if(lp == std::string::npos) {
        name = temp_name;
    } else {
        name = temp_name.substr(leader.size());
    }

    std::replace(name.begin(), name.end(), '/', '_');

    metadata = get_model_metadata_from_name(name);
}

void PoreModel::write(const std::string filename, const std::string modelname) const
{
    std::string outmodelname = modelname;
    if(modelname.empty())
        outmodelname = name;

    std::ofstream writer(filename);
    writer << "#model_name\t" << outmodelname << std::endl;
    writer << "#type\t" << this->type << std::endl;
    writer << "#kit\t" << this->metadata.get_kit_name() << std::endl;
    writer << "#strand\t" << this->metadata.get_strand_model_name() << std::endl;
    writer << "#shift_offset\t" << this->shift_offset << std::endl;
    writer << "#scale_offset\t" << this->scale_offset << std::endl;

    std::string curr_kmer(k, this->pmalphabet->base(0));
    for(size_t ki = 0; ki < this->states.size(); ++ki) {
        writer << curr_kmer << "\t" << this->states[ki].level_mean << "\t" << this->states[ki].level_stdv << "\t"
               << this->states[ki].sd_mean << "\t" << this->states[ki].sd_stdv << std::endl;
        this->pmalphabet->lexicographic_next(curr_kmer);
    }
    writer.close();
}

void PoreModel::update_states( const PoreModel &other )
{
    k = other.k;
    pmalphabet = other.pmalphabet;
    shift += other.shift_offset;
    scale += other.scale_offset;
    update_states( other.states );
}

void PoreModel::update_states( const std::vector<PoreModelStateParams> &otherstates )
{
    states = otherstates;
    if (is_scaled) {
        bake_gaussian_parameters();
    }
}

void PoreModel::set_metadata(const std::string& kit, const std::string& strand)
{
    assert(!kit.empty());
    assert(!strand.empty());

    if(kit == "SQK006") {
        this->metadata.kit = KV_SQK006;
    } else if(kit == "r9_250bps") {
        this->metadata.kit = KV_R9_250BPS;
    } else if(kit == "r9.4_450bps") {
        this->metadata.kit = KV_R9_4_450BPS;
    } else {
        fprintf(stderr, "Error, unrecognized model kit %s\n", kit.c_str());
        exit(EXIT_FAILURE);
    }

    if(strand == "template") {
        this->metadata.model_idx = 0;
        this->metadata.strand_idx = 0;
    } else if(strand == "complement.pop1") {
        this->metadata.model_idx = 1;
        this->metadata.strand_idx = 1;
    } else if(strand == "complement.pop2") {
        this->metadata.model_idx = 2;
        this->metadata.strand_idx = 1;
    } else {
        fprintf(stderr, "Error, unrecognized model strand %s\n", strand.c_str());
        exit(EXIT_FAILURE);
    }
}
