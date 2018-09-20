//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_variants -- tools for calling variants
//
#include <algorithm>
#include <map>
#include <iterator>
#include <iomanip>
#include "nanopolish_profile_hmm.h"
#include "nanopolish_variant.h"
#include "nanopolish_haplotype.h"
#include "nanopolish_model_names.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_variant_db.h"

//#define DEBUG_HAPLOTYPE_SELECTION 1

std::string Variant::make_vcf_header_key_value(const std::string& key, const std::string& value)
{
    std::stringstream ss;
    ss << "##" << key << "=" << value;
    return ss.str();
}

std::string Variant::make_vcf_tag_string(const std::string& tag,
                                         const std::string& id,
                                         int count,
                                         const std::string& type,
                                         const std::string& description)
{
    std::stringstream ss;
    ss << "##" << tag << "=<ID=" << id << ",Number=" << count << ",Type="
       << type << ",Description=\"" << description << "\">";
    return ss.str();
}

void Variant::write_vcf_header(FILE* fp,
                               const std::vector<std::string>& header_lines)
{

    fprintf(fp, "##fileformat=VCFv4.2\n");
    for(const std::string& line : header_lines) {
        fprintf(fp, "%s\n", line.c_str());
    }
    fprintf(fp, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample\n");
}

// return a new copy of the string with gap symbols removed
std::string remove_gaps(const std::string& str)
{
    std::string ret = str;
    ret.erase( std::remove(ret.begin(), ret.end(), '-'), ret.end());
    return ret;
}

std::vector<Variant> read_variants_from_file(const std::string& filename)
{
    std::vector<Variant> out;
    std::ifstream infile(filename);
    std::string line;
    while(getline(infile, line)) {

        // skip headers
        if(line[0] == '#') {
            continue;
        }

        // parse variant
        Variant v(line);
        v.info.clear();
        out.push_back(v);
    }
    return out;
}

std::vector<Variant> read_variants_for_region(const std::string& filename,
                                              const std::string& contig,
                                              int region_start,
                                              int region_end)
{
    std::vector<Variant> out;
    std::ifstream infile(filename);
    std::string line;
    while(getline(infile, line)) {

        // skip headers
        if(line[0] == '#') {
            continue;
        }

        // parse variant
        Variant v(line);

        if(v.ref_name == contig &&
           (int)v.ref_position >= region_start &&
           (int)v.ref_position <= region_end)
        {
            v.info.clear();
            out.push_back(v);
        }
    }
    return out;
}

void filter_variants_by_count(std::vector<Variant>& variants, int min_count)
{
    std::map<std::string, std::pair<Variant, int>> map;

    for(size_t i = 0; i < variants.size(); ++i) {
        std::string key = variants[i].key();
        auto iter = map.find(key);
        if(iter == map.end()) {
            map.insert(std::make_pair(key, std::make_pair(variants[i], 1)));
        } else {
            iter->second.second += 1;
        }
    }

    variants.clear();

    for(auto iter = map.begin(); iter != map.end(); ++iter) {
        Variant& v = iter->second.first;
        if(iter->second.second >= min_count) {
            v.add_info("BaseCalledReadsWithVariant", iter->second.second);
            variants.push_back(v);
        }
    }
}

void filter_out_non_snp_variants(std::vector<Variant>& variants)
{
    std::vector<Variant> tmp;
    for(size_t i = 0; i < variants.size(); ++i) {
        bool is_snp = variants[i].ref_seq.length() == 
                      variants[i].alt_seq.length();
        if(is_snp) {
            tmp.push_back(variants[i]);
        }
    }
    variants.swap(tmp);
}

// this doesn't handle triallele ("1/2") genotypes, yet
std::string make_genotype(size_t alt_alleles, size_t ploidy)
{
    assert(alt_alleles <= ploidy);
    std::string out(ploidy + ploidy - 1, '/');
    int num_ref = ploidy - alt_alleles;
    for(size_t pos = 0; pos < out.size(); pos += 2) {
        out[pos] = (pos / 2) < num_ref ? '0' : '1';
    }
    return out;
}

//
std::vector<HMMInputSequence> generate_methylated_alternatives(const HMMInputSequence& sequence, 
                                                               const std::vector<std::string>& methylation_types)
{
    // Make methylated versions of each input sequence
    std::vector<HMMInputSequence> out;
    out.push_back(sequence);

    for(size_t i = 0; i < methylation_types.size(); ++i) {
        const Alphabet* alphabet = get_alphabet_by_name(methylation_types[i]);
        std::string methylated = alphabet->methylate(sequence.get_sequence());

        // Is there a methylated version?
        if(methylated != sequence.get_sequence()) {
            out.emplace_back(methylated, alphabet);
        }
    }
    return out;
}


//
void score_variant_group(VariantGroup& variant_group,
                         Haplotype base_haplotype, 
                         const std::vector<HMMInputData>& input,
                         const int max_haplotypes,
                         const int ploidy,
                         const bool genotype_all_input_variants,
                         const uint32_t alignment_flags,
                         const std::vector<std::string>& methylation_types)
{
    size_t num_variants = variant_group.get_num_variants();

    // Determine the maximum number of variants we can jointly test
    // without exceeding the maximum number of haplotypes
    size_t sum_num_haplotypes = 0;
    size_t max_r = 1;
    
    while(max_r <= num_variants) {
        size_t num_haplotypes_r = nChoosek(num_variants, max_r);
        if(num_haplotypes_r + sum_num_haplotypes < max_haplotypes) {
            sum_num_haplotypes += num_haplotypes_r;
        } else {
            break;
        }
        max_r += 1;
    }
    max_r -= 1;

    if(max_r != num_variants) {
        fprintf(stderr, "Number of variants in span (%lu) would exceed max-haplotypes. Variants may be missed. Consider running with a higher value of max-haplotypes!\n", num_variants);
    }

    // Construct haplotypes (including the base haplotype with no variants)
    // Track the variant combination ID within the group
    std::vector<std::pair<Haplotype, size_t>> haplotypes;

    for(size_t r = 0; r <= max_r; ++r) {

        Combinations combinations(num_variants, r);
        while(!combinations.done()) {
            VariantCombination vc(combinations.get());

            // Apply variants to haplotype
            Haplotype current_haplotype = base_haplotype;
            bool good_haplotype = current_haplotype.apply_variants(variant_group.get_variants(vc));
            if(good_haplotype) {
                size_t vc_idx = variant_group.add_combination(vc);
                haplotypes.push_back(std::make_pair(current_haplotype, vc_idx));
            }
            combinations.next();
        }
    }

    std::vector<std::string> read_ids;
    for(size_t i = 0; i < input.size(); ++i) {
        std::stringstream ss;
        ss << input[i].read->read_name << ":" << input[i].strand;
        read_ids.push_back(ss.str());
    }
  
    #pragma omp parallel for
    for(size_t ri = 0; ri < input.size(); ++ri) {
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            const auto& current = haplotypes[hi];

            // Expand the haplotype to contain all representations of this sequence by adding methylation
            std::vector<HMMInputSequence> sequences = generate_methylated_alternatives(current.first.get_sequence(), methylation_types);
            double score = profile_hmm_score_set(sequences, input[ri], alignment_flags);
            
            #pragma omp critical
            {
                variant_group.set_combination_read_score(current.second, read_ids[ri], score);
//                set(read_haplotype_scores, ri, hi, score);
//                read_sum[ri] = add_logs(read_sum[ri], score);
            }
        }
    }

#if 0
#endif
}

std::vector<Variant> simple_call(VariantGroup& variant_group,
                                 const int ploidy,
                                 const bool genotype_all_input_variants)
{
    // Select the best set of haplotypes that maximizes the probability of the data
    double base_score = -INFINITY;
    double best_score = -INFINITY;
    std::vector<size_t> best_set;
    std::vector<size_t> base_set;

    double log_2 = log(2);

#ifdef DEBUG_HAPLOTYPE_SELECTION 
    fprintf(stderr, "Selecting haplotypes\n");
#endif
    
    // Get read data for this group
    const std::vector< std::pair<std::string, double>> group_reads = variant_group.get_read_sum_scores();

    // Skip groups that only have one possibility (these are typically malformed VCF records)
    size_t variant_combos_in_group = variant_group.get_num_combinations();
    if (variant_combos_in_group <= 1) {
        return std::vector<Variant>();
    }

    Combinations vc_sets(variant_combos_in_group, ploidy, CO_WITH_REPLACEMENT);
    while(!vc_sets.done()) {

        // The current combination is represented as a vector of haplotype IDs
        std::vector<size_t> current_set = vc_sets.get();
        
        // Check if the current set consists of entirely of haplotypes without variants
        bool is_base_set = true;
        for(size_t i = 0; i < current_set.size(); ++i) {
            const VariantCombination& vc = variant_group.get_combination(current_set[i]);
            is_base_set = is_base_set && (variant_group.get_variants(vc).size() == 0);
        }

        double set_score = 0.0f;
        std::vector<double> read_support(current_set.size(), 0.0f);

        for(size_t ri = 0; ri < group_reads.size(); ++ri) {
            const std::string& read_id = group_reads[ri].first;

            double set_sum = -INFINITY;
            for(size_t j = 0; j < current_set.size(); ++j) {
                size_t vc_id = current_set[j];
                double rhs = variant_group.get_combination_read_score(vc_id, read_id);
                /*
                fprintf(stderr, "\t\tread-haplotype: %s %zu %s %.2lf\n", read_id.c_str(), 
                                                                         variant_group.get(0).ref_position, // hack
                                                                         variant_group.get_vc_allele_string(vc_id).c_str(), rhs); 
                */
                set_sum = add_logs(set_sum, rhs - log_2);
                read_support[j] += exp(rhs - group_reads[ri].second);
            }

            /*
            fprintf(stderr, "\t\tread-genotype: %s %zu %.2lf\n", read_id.c_str(), 
                                                                         variant_group.get(0).ref_position, // hack
                                                                         set_sum); 
            */
            set_score += set_sum;
        }
        
        if(is_base_set) {
            base_score = set_score;
            base_set = current_set;
        }
        
#ifdef DEBUG_HAPLOTYPE_SELECTION 
        fprintf(stderr, "Current set score: %.5lf\t", set_score);
        for(size_t i = 0; i < current_set.size(); ++i) {
            fprintf(stderr, "\t%zu:%.2lf", current_set[i], read_support[i]);
        }
        fprintf(stderr, "\n");
#endif

        if(set_score > best_score) {
            best_score = set_score;
            best_set = current_set;
        }
        vc_sets.next();
    }

    // TODO: set an appropriate threshold
    if(best_score - base_score < 20) {
        best_set = base_set;
    }

    // Calculate the number of reads that support each variant allele
    std::vector<double> read_variant_support(variant_group.get_num_variants(), 0.0f);
    for(size_t vc_id = 0; vc_id < variant_group.get_num_combinations(); ++vc_id) {

        const VariantCombination& vc = variant_group.get_combination(vc_id);
        
        for(size_t ri = 0; ri < group_reads.size(); ++ri) {
            const std::string& read_id = group_reads[ri].first;
            double read_sum = group_reads[ri].second;
            double read_haplotype_score = variant_group.get_combination_read_score(vc_id, read_id);
            double posterior_read_from_haplotype = exp(read_haplotype_score - read_sum);

            for(size_t var_idx = 0; var_idx < vc.get_num_variants(); ++var_idx) {
                read_variant_support[vc.get_variant_id(var_idx)] += posterior_read_from_haplotype;
            }
        }
    }

    std::vector<Variant> output_variants;
    for(size_t vi = 0; vi < variant_group.get_num_variants(); vi++) {

        size_t var_count = 0;
        for(size_t j = 0; j < best_set.size(); ++j) {
            const VariantCombination& curr_vc = variant_group.get_combination(best_set[j]);

            for(size_t k = 0; k < curr_vc.get_num_variants(); ++k) {
                var_count += curr_vc.get_variant_id(k) == vi;
            }
        }
        
        if( !(genotype_all_input_variants || var_count > 0)) {
            continue;
        }

        Variant v = variant_group.get(vi);
        if(var_count > 0) {
            v.quality = best_score - base_score;
        } else {
            v.quality = 0.0;
        }
        v.add_info("TotalReads", group_reads.size());
        v.add_info("AlleleCount", var_count);
        v.add_info("SupportFraction", read_variant_support[vi] / group_reads.size());
        v.genotype = make_genotype(var_count, ploidy);
        output_variants.push_back(v);
    }

    return output_variants;

}

std::string prettyprint_multi(std::vector<size_t>& group_permutation_indices, 
                              const std::vector<SizeTVecVec>& variant_group_permutations, 
                              std::vector<const VariantGroup*> all_groups)
{
    std::stringstream ss;
    for(size_t group_idx = 0; group_idx < group_permutation_indices.size(); ++group_idx) {
        size_t group_perm_idx = group_permutation_indices[group_idx];
        const auto& vc_vec = variant_group_permutations[group_idx][group_perm_idx]; 
        for(size_t vc_idx = 0; vc_idx < vc_vec.size(); ++vc_idx) {
            size_t id = vc_vec[vc_idx];
            char sep = vc_idx == vc_vec.size() - 1 ? ' ' : '/';
            ss << all_groups[group_idx]->get_vc_allele_string(id) << sep;
        }
    }
    std::string out = ss.str();
    return out.substr(0, out.size() - 1);
}

std::string prettyprint_haplotype(const std::vector<size_t>& haplotype, 
                                  std::vector<const VariantGroup*> all_groups)
{
    std::stringstream ss;
    for(size_t group_idx = 0; group_idx < haplotype.size(); ++group_idx) {
        const auto& vc_idx = haplotype[group_idx];
        const char* sep = group_idx == haplotype.size() - 1 ? "" : ",";
        ss << all_groups[group_idx]->get_vc_allele_string(vc_idx) << sep;
    }
    return ss.str();
}

std::string prettyprint_genotype(std::vector<size_t>& genotype, 
                                 const SizeTVecVec& haplotypes, 
                                 std::vector<const VariantGroup*> all_groups)
{
    std::stringstream ss;
    for(size_t gt_idx = 0; gt_idx < genotype.size(); ++gt_idx) {
        size_t hap_idx = genotype[gt_idx];

        // The haplotype should have one variant combo per group
        assert(haplotypes[hap_idx].size() == all_groups.size());

        for(size_t group_idx = 0; group_idx < haplotypes[hap_idx].size(); ++group_idx) {
            const auto& vc_idx = haplotypes[hap_idx][group_idx];
            const char* sep = group_idx == haplotypes[hap_idx].size() - 1 ? "" : ",";
            ss << all_groups[group_idx]->get_vc_allele_string(vc_idx) << sep;
        }
        ss << "/";
    }
    std::string out = ss.str();
    return out.substr(0, out.size() - 1);
}

typedef std::vector<size_t> Genotype;

std::vector<Variant> multi_call(VariantGroup& variant_group,
                                std::vector<const VariantGroup*> neighbor_groups,
                                const int ploidy,
                                const bool genotype_all_input_variants)
{
    // Select the best set of haplotypes that maximizes the probability of the data
    double base_score = -INFINITY;
    double best_score = -INFINITY;
    std::vector<size_t> best_set;
    std::vector<size_t> base_set;

#ifdef DEBUG_HAPLOTYPE_SELECTION 
    fprintf(stderr, "Selecting haplotypes\n");
#endif

    std::vector<const VariantGroup*> all_groups;
    all_groups.push_back(&variant_group);
    all_groups.insert(all_groups.end(), neighbor_groups.begin(), neighbor_groups.end());
    
    // Build an index of each variant combination for each variant group
    SizeTVecVec variant_combinations_by_group(all_groups.size());
    for(size_t gi = 0; gi < all_groups.size(); ++gi) {
        std::vector<size_t> vc_ids;
        for(size_t vci = 0; vci < all_groups[gi]->get_num_combinations(); ++vci) {
            variant_combinations_by_group[gi].push_back(vci);
        }
    }

    // Build haplotypes by generating all permutations of the variant combos for each group
    SizeTVecVec haplotypes = cartesian_product(variant_combinations_by_group);
    
    // get the read data for the variant group that we are genotyping
    const std::vector< std::pair<std::string, double>> group_reads = variant_group.get_read_sum_scores();

    // Score each haplotype
    DoubleMatrix read_haplotype_scores;
    allocate_matrix(read_haplotype_scores, group_reads.size(), haplotypes.size());
    
    // Calculate and store read-haplotype scores
    for(size_t ri = 0; ri < group_reads.size(); ++ri) {
        auto read_id = group_reads[ri].first;
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {

            const auto& haplotype = haplotypes[hi];

            // The haplotype should have one variant combo per group
            assert(haplotype.size() == all_groups.size());
            double hap_sum = 0.0f;
            for(size_t group_idx = 0; group_idx < haplotype.size(); ++group_idx) {
                const auto& vc_idx = haplotype[group_idx];
                hap_sum += all_groups[group_idx]->get_combination_read_score(vc_idx, read_id);
            }

            set(read_haplotype_scores, ri, hi, hap_sum);
        }
    }

    

    // Dindel EM model
    // Calculate expectation of read-haplotype indicator variables
    DoubleMatrix z;
    allocate_matrix(z, group_reads.size(), haplotypes.size());
    for(size_t ri = 0; ri < group_reads.size(); ++ri) {
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            set(z, ri, hi, 0.5); // doEM initializes to 0.5, should be 1/haplotypes.size()?
        }
    }
    
    // log of haplotype frequencies
    std::vector<double> pi(haplotypes.size(), log(1.0 / haplotypes.size()));

    // counts of each haplotype
    std::vector<double> nk(haplotypes.size(), 0.0f);
    
    // run EM
    size_t iterations = 0;
    while(1) {

        for(size_t i = 0; i < haplotypes.size(); ++i) {
            nk[i] = 0.0;
        }

        for(size_t ri = 0; ri < group_reads.size(); ++ri) {

            // responsibility
            double lognorm = -INFINITY;
            for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
                set(z, ri, hi, pi[hi] + get(read_haplotype_scores, ri, hi));
                lognorm = add_logs(lognorm, get(z, ri, hi));
            }

            // normalize
            for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
                double t = get(z, ri, hi);
                double et = exp(t - lognorm);
                set(z, ri, hi, et);
                nk[hi] += et;
            }
        }

        // frequencies
        double zh = 0.0f;
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            zh += nk[hi];
        }

        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            pi[hi] = log(nk[hi] / zh);
        }

        // debug output
        for(size_t ri = 0; ri < group_reads.size(); ++ri) {
            fprintf(stderr, "read-haplotype indicator - %s\t", group_reads[ri].first.c_str());
            for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
                std::string hap_str = prettyprint_haplotype(haplotypes[hi], all_groups);
                fprintf(stderr, "%s: %.3lf ", hap_str.c_str(), get(z, ri, hi));
            }
            fprintf(stderr, "\n");
        }

        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            std::string hap_str = prettyprint_haplotype(haplotypes[hi], all_groups);
            fprintf(stderr, "hap[%zu]: %s read count: %.2lf\n", hi, hap_str.c_str(), nk[hi]);
        }

        fprintf(stderr, "haplotype frequencies:");
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            std::string hap_str = prettyprint_haplotype(haplotypes[hi], all_groups);
            fprintf(stderr, "%s:%.3lf ", hap_str.c_str(), exp(pi[hi]));
        }
        fprintf(stderr, "\n");
        
        if(iterations++ > 2) {
            break;
        }
    }

    // Build genotypes by selecting ploidy haplotypes (with replacement)
    Combinations genotype_combos(haplotypes.size(), ploidy, CO_WITH_REPLACEMENT);
    SizeTVecVec genotypes;
    while(!genotype_combos.done()) {
        genotypes.push_back(genotype_combos.get());
        genotype_combos.next();
    }
   
    // Score genotypes
    double log_2 = log(2.0f);
    std::vector<double> scores(genotypes.size(), 0.0);
    for(size_t i = 0; i < genotypes.size(); ++i) {
        const auto& genotype = genotypes[i];

        // Score all reads against this genotype
        for(size_t ri = 0; ri < group_reads.size(); ++ri) {
            
            double read_sum = -INFINITY;

            for(size_t gt_idx = 0; gt_idx < genotype.size(); ++gt_idx) {
                size_t hi = genotype[gt_idx];
                double read_hap_score = get(read_haplotype_scores, ri, hi);
                const auto& haplotype = haplotypes[genotype[gt_idx]];
                std::string hap_str = prettyprint_haplotype(haplotype, all_groups);
                fprintf(stderr, "\t\t%s %s %.2lf\n", group_reads[ri].first.c_str(), hap_str.c_str(), read_hap_score);
                read_sum = add_logs(read_sum, read_hap_score - log_2);
            }
            scores[i] += read_sum;
        }
    }

    for(size_t si = 0; si < scores.size(); ++si) {
        std::string perm_pp = prettyprint_genotype(genotypes[si], haplotypes, all_groups);
        fprintf(stderr, "perm[%zu]: %s %.2lf\n", si, perm_pp.c_str(), scores[si]);
    }

#if 0
    // TODO: set an appropriate threshold
    if(best_score - base_score < 5) {
        best_set = base_set;
    }

    std::vector<Variant> output_variants;
    for(size_t vi = 0; vi < variant_group.get_num_variants(); vi++) {

        size_t var_count = 0;
        for(size_t j = 0; j < best_set.size(); ++j) {
            const VariantCombination& curr_vc = variant_group.get_combination(best_set[j]);

            for(size_t k = 0; k < curr_vc.get_num_variants(); ++k) {
                var_count += curr_vc.get_variant_id(k) == vi;
            }
        }
        
        if( !(genotype_all_input_variants || var_count > 0)) {
            continue;
        }

        Variant v = variant_group.get(vi);
        if(var_count > 0) {
            v.quality = best_score - base_score;
        } else {
            v.quality = 0.0;
        }
        v.add_info("TotalReads", group_reads.size());
        v.add_info("AlleleCount", var_count);
        v.genotype = make_genotype(var_count, ploidy);
        output_variants.push_back(v);
    }
#endif
    std::vector<Variant> output_variants;

    free_matrix(read_haplotype_scores);
    free_matrix(z);
    return output_variants;
}

//
Variant score_variant_thresholded(const Variant& input_variant,
                                  Haplotype base_haplotype, 
                                  const std::vector<HMMInputData>& input,
                                  const uint32_t alignment_flags,
                                  const uint32_t score_threshold,
                                  const std::vector<std::string>& methylation_types)
{

    Variant out_variant = input_variant;

    Haplotype variant_haplotype = base_haplotype;
    variant_haplotype.apply_variant(input_variant);

    // Make methylated versions of each input sequence
    std::vector<HMMInputSequence> base_sequences = generate_methylated_alternatives(base_haplotype.get_sequence(), methylation_types);
    std::vector<HMMInputSequence> variant_sequences = generate_methylated_alternatives(variant_haplotype.get_sequence(), methylation_types);
    
    double total_score = 0.0f;
    #pragma omp parallel for
    for(size_t j = 0; j < input.size(); ++j) {

        if(fabs(total_score) < score_threshold) {

            // Calculate scores using the base nucleotide model
            double base_score = profile_hmm_score_set(base_sequences, input[j], alignment_flags);
            double variant_score = profile_hmm_score_set(variant_sequences, input[j], alignment_flags);

            #pragma omp atomic
            total_score += (variant_score - base_score);
        }
    }

    out_variant.quality = total_score;
    return out_variant;
}

void annotate_variants_with_all_support(std::vector<Variant>& input, const AlignmentDB& alignments, int min_flanking_sequence, const uint32_t alignment_flags)
{
    Haplotype ref_haplotype(alignments.get_region_contig(), alignments.get_region_start(), alignments.get_reference());

    for(size_t vi = 0; vi < input.size(); vi++) {
        Variant& v = input[vi];

        std::string ref_name = v.ref_name;
        int calling_start = v.ref_position - min_flanking_sequence;
        int calling_end = v.ref_position + min_flanking_sequence;

        // Construct haplotype set for A,C,G,T at this position
        std::vector<Haplotype> haplotypes;
        Haplotype calling_haplotype = ref_haplotype.substr_by_reference(calling_start, calling_end);

        Variant tmp_variant = v;
        for(size_t bi = 0; bi < 4; bi++) {
            Haplotype variant_haplotype = calling_haplotype;
            tmp_variant.alt_seq = "ACGT"[bi];
            variant_haplotype.apply_variant(tmp_variant);
            haplotypes.push_back(variant_haplotype);
        }

        // Output
        std::vector<double> support_fraction(4);

        std::vector<HMMInputData> event_sequences =
            alignments.get_event_subsequences(ref_name, calling_start, calling_end);

        for(size_t ri = 0; ri < event_sequences.size(); ++ri) {
            double sum_score = 0.0;
            std::vector<double> scores;
            for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
                double s = profile_hmm_score(haplotypes[hi].get_sequence(), event_sequences[ri], alignment_flags);
                scores.push_back(s);
                if(hi > 0) {
                    sum_score = add_logs(s, sum_score);
                } else {
                    sum_score = s;
                }
            }

            for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
                support_fraction[hi] += exp(scores[hi] - sum_score);
            }
        }

        std::stringstream ss;
        ss << std::fixed << std::setprecision(3);
        for(size_t hi = 0; hi < haplotypes.size(); ++hi) {
            support_fraction[hi] /= event_sequences.size();
            ss << support_fraction[hi];
            if(hi != haplotypes.size() - 1) {
                ss << ",";
            }
        }
        v.add_info("SupportFractionByBase", ss.str());
    }
}
