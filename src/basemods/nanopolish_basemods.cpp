//---------------------------------------------------------
// Copyright 2021 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_basemods -- utility functions for working with
// base modification data
//
#include "nanopolish_basemods.h"
#include <map>
#include <string>
#include "nanopolish_alignment_db.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_bam_utils.h"


std::string retrieve_unambiguous_reference_sequence(const faidx_t* fai, const bam_hdr_t* hdr, const bam1_t* record)
{
    std::string contig = hdr->target_name[record->core.tid];
    int ref_start_pos = record->core.pos;
    int ref_end_pos =  bam_endpos(record);

    // Extract the reference sequence for this region
    int fetched_len = 0;
    assert(ref_end_pos >= ref_start_pos);
    std::string ref_seq = get_reference_region_ts(fai, contig.c_str(), ref_start_pos,
                                                  ref_end_pos, &fetched_len);

    // Remove non-ACGT bases from this reference segment
    ref_seq = gDNAAlphabet.disambiguate(ref_seq);
    return ref_seq;
}

void get_modification_symbols(const Alphabet* alphabet, char& unmodified_symbol, char modified_symbol)
{
    // determine the symbol for the unmodified and modified base
    assert(alphabet->num_recognition_sites() == 1);
    unmodified_symbol = 'N';
    modified_symbol = METHYLATED_SYMBOL;
    const char* modified_motif = alphabet->get_recognition_site_methylated(0);
    for(size_t i = 0; i < alphabet->recognition_length(); ++i) {
        if(modified_motif[i] == modified_symbol) {
            assert(unmodified_symbol == 'N');
            unmodified_symbol = alphabet->get_recognition_site(0)[i];
        }
    }
    assert(unmodified_symbol != 'N');
}

void calculate_call_vectors(const std::map<int, ScoredSite>& calls, 
                            const Alphabet* alphabet, 
                            std::vector<size_t>& call_reference_positions,
                            std::vector<uint8_t>& call_probabilities)
{
    for(auto iter : calls) {
        const ScoredSite& call = iter.second;
        const std::string& seq = call.sequence;

        // convert the positions we've called at to Ms
        std::string m_seq = alphabet->methylate(seq);

        // start_position gives the position of the first site in the group
        // we need to record an offset here so we don't count the flank
        int start_pos = call.start_position;
        size_t flank_offset = m_seq.find_first_of(METHYLATED_SYMBOL);
        assert(flank_offset != std::string::npos);
        
        // this is shared across all CGs in the site
        double methylation_probability = exp(call.ll_methylated[0]) / ( exp(call.ll_methylated[0]) + exp(call.ll_unmethylated[0]) );
        uint8_t methylation_probability_code = std::min(255, (int)(methylation_probability * 255));

        for(size_t j = 0; j < m_seq.size(); j++) {
            if(m_seq[j] == METHYLATED_SYMBOL) {
                int reference_position = start_pos + j - flank_offset;
                call_reference_positions.push_back(reference_position);
                call_probabilities.push_back(methylation_probability_code);
            }
        }
    }
}

std::string generate_mm_tag(char unmodified_symbol, const std::string& sequence, const std::vector<size_t>& call_seq_indices)
{
    std::string delta_str;
    delta_str += unmodified_symbol;
    delta_str += "+m";

    int count_start = 0;
    for(size_t call_index = 0; call_index < call_seq_indices.size(); ++call_index) {
        // Count the number of unmodified bases preceding this position
        int count = 0;
        for(size_t j = count_start; j < call_seq_indices[call_index]; ++j) {
            count += sequence[j] == unmodified_symbol;
        }
        
        // this intentionally adds a delimiter to the first element
        delta_str += ',';
        delta_str += std::to_string(count);

        count_start = call_seq_indices[call_index] + 1;
    }

    delta_str += ';';
    return delta_str;
}

bam1_t* create_modbam_record(const bam1_t* record,
                             const std::map<int, ScoredSite>& calls,
                             const MethylationCallingParameters& calling_parameters)
{
    
    char unmodified_symbol = 'N';
    char modified_symbol = 'N';
    get_modification_symbols(calling_parameters.alphabet, unmodified_symbol, modified_symbol);

    // Convert the calls into a vector of reference positions and probabilities
    std::vector<size_t> call_reference_positions; 
    std::vector<uint8_t> call_reference_probabilities; 

    calculate_call_vectors(calls, calling_parameters.alphabet, call_reference_positions, call_reference_probabilities);

    // create a map from reference position -> read position from the CIGAR string in the BAM
    bool rc = bam_is_rev(record);
    SequenceAlignmentRecord alignment(record);

    std::string original_sequence = !rc ? alignment.sequence : gDNAAlphabet.reverse_complement(alignment.sequence);

    std::map<size_t, size_t> reference_to_read_map;
    for(auto& aligned_pair : alignment.aligned_bases) {
        reference_to_read_map[aligned_pair.ref_pos] = !rc ? aligned_pair.read_pos : original_sequence.length() - aligned_pair.read_pos - 1;
    }
    
    // TODO: support/test other modes
    int strand_offset = !rc ? 0 : 1;
    assert(calling_parameters.alphabet->get_name() == "cpg");

    // convert the call reference positions into a vector of indices into the SEQ record of the BAM (which is held in SequenceAlignmentRecord)
    // NB some calls may get dropped so we need to make a new probabilities vector too
    std::vector<size_t> call_seq_indices;
    std::vector<uint8_t> call_seq_probabilities;
    for(size_t i = 0; i < call_reference_positions.size(); ++i) {

        // the rc offset is to account for that fact that when sequencing a CG on the opposite strand
        // the read base we're interested in is aligned to the G

        // ref       ATTACGTTTA
        // rc(read)  ATTACGTTTA
        //                ^
        size_t reference_position = call_reference_positions[i] + strand_offset;
        auto iter = reference_to_read_map.find(reference_position);
        if(iter != reference_to_read_map.end()) {
            size_t read_index = iter->second;
            assert(read_index < original_sequence.size());
            if(original_sequence[read_index] == unmodified_symbol) {
                call_seq_indices.push_back(read_index);
                call_seq_probabilities.push_back(call_reference_probabilities[i]);
            }
        }
    }

    // reverse arrays when SEQ is reverse complemented so we proceed in the original direction
    if(rc) {
        std::reverse(call_seq_indices.begin(), call_seq_indices.end());
        std::reverse(call_seq_probabilities.begin(), call_seq_probabilities.end());
    }

    std::string delta_str = generate_mm_tag(unmodified_symbol, original_sequence, call_seq_indices);

    bam1_t* mod_record = bam_dup1(record);

    int status = bam_aux_update_str(mod_record, "Mm", delta_str.size() + 1, delta_str.c_str());
    assert(status == 0);
    
    status = bam_aux_update_array(mod_record, "Ml", 'C', call_seq_probabilities.size(), call_seq_probabilities.data());
    assert(status == 0);
    return mod_record;
}

bam1_t* create_reference_modbam_record(const faidx_t* fai,
                                       const bam_hdr_t* hdr,
                                       const bam1_t* record,
                                       const std::map<int, ScoredSite>& calls,
                                       const MethylationCallingParameters& calling_parameters)
{
    char unmodified_symbol = 'N';
    char modified_symbol = 'N';
    get_modification_symbols(calling_parameters.alphabet, unmodified_symbol, modified_symbol);

    std::vector<size_t> call_reference_positions; 
    std::vector<uint8_t> call_reference_probabilities; 
    calculate_call_vectors(calls, calling_parameters.alphabet, call_reference_positions, call_reference_probabilities);

    // convert the call reference positions into a vector of indices into ref_seq
    std::vector<size_t> indices; 

    int ref_start_pos = record->core.pos;
    for(size_t i = 0; i < call_reference_positions.size(); ++i) {
        indices.push_back(call_reference_positions[i] - ref_start_pos);
    }

    std::string ref_seq = retrieve_unambiguous_reference_sequence(fai, hdr, record);
    std::string delta_str = generate_mm_tag(unmodified_symbol, ref_seq, indices);

    // rewrite SEQ/cigar with reference sequence
    bam1_t* out_record = bam_init1();

    // basic stats
    out_record->core.tid = record->core.tid;
    out_record->core.pos = record->core.pos;
    out_record->core.qual = record->core.qual;
    out_record->core.flag = 0;
    out_record->core.bin = record->core.bin;

    // no read pairs
    out_record->core.mtid = -1;
    out_record->core.mpos = -1;
    out_record->core.isize = 0;

    std::string read_name = bam_get_qname(record);
    std::string outqual(ref_seq.length(), 30);

    std::vector<uint32_t> cigar;
    uint32_t cigar_op = ref_seq.size() << BAM_CIGAR_SHIFT | BAM_CMATCH;
    cigar.push_back(cigar_op);
    write_bam_vardata(out_record, read_name, cigar, ref_seq, outqual);
    
    fprintf(stderr, "%s STR %s\n", bam_get_qname(record), delta_str.c_str());
    int status = bam_aux_update_str(out_record, "Mm", delta_str.size() + 1, delta_str.c_str());
    assert(status == 0);
    
    status = bam_aux_update_array(out_record, "Ml", 'C', call_reference_probabilities.size(), call_reference_probabilities.data());
    assert(status == 0);

    return out_record;
}

// Test motif sites in this read for methylation
void calculate_methylation_for_read(const OutputHandles& handles,
                                    MethylationCallingResult& result,
                                    SquiggleRead& sr,
                                    const MethylationCallingParameters& calling_parameters,
                                    const faidx_t* fai,
                                    const bam_hdr_t* hdr,
                                    const bam1_t* record,
                                    size_t read_idx,
                                    int region_start,
                                    int region_end)
{

    // Insert this bam record into the result structure to initialize a new site score map
    assert(record != NULL);

    std::map<int, ScoredSite>* site_score_map = NULL;
    #pragma omp critical
    {
        site_score_map = &result[record];
    }
    
    // Retrieve reference
    std::string contig = hdr->target_name[record->core.tid];
    int ref_start_pos = record->core.pos;
    int ref_end_pos =  bam_endpos(record);

    // Extract the reference sequence for this region
    int fetched_len = 0;
    assert(ref_end_pos >= ref_start_pos);
    std::string ref_seq = get_reference_region_ts(fai, contig.c_str(), ref_start_pos,
                                                  ref_end_pos, &fetched_len);

    // Remove non-ACGT bases from this reference segment
    ref_seq = gDNAAlphabet.disambiguate(ref_seq);

    for(size_t strand_idx = 0; strand_idx < NUM_STRANDS; ++strand_idx) {
        if(!sr.has_events_for_strand(strand_idx)) {
            continue;
        }

        size_t k = sr.get_model_k(strand_idx);

        // check if there is a motif model for this strand
        if(!PoreModelSet::has_model(sr.get_model_kit_name(strand_idx),
                                    calling_parameters.methylation_type,
                                    sr.get_model_strand_name(strand_idx),
                                    k))
        {
            continue;
        }

        // Build the event-to-reference map for this read from the bam record
        SequenceAlignmentRecord seq_align_record(record);
        EventAlignmentRecord event_align_record(&sr, strand_idx, seq_align_record);

        std::vector<double> site_scores;
        std::vector<int> site_starts;
        std::vector<int> site_ends;
        std::vector<int> site_count;

        // Scan the sequence for motifs
        std::vector<int> motif_sites;
        assert(ref_seq.size() != 0);
        for(size_t i = 0; i < ref_seq.size() - 1; ++i) {
            if(calling_parameters.alphabet->is_motif_match(ref_seq, i))
                motif_sites.push_back(i);
        }

        // Batch the motifs together into groups that are separated by some minimum distance
        std::vector<std::pair<int, int>> groups;

        size_t curr_idx = 0;
        while(curr_idx < motif_sites.size()) {
            // Find the endpoint of this group of sites
            size_t end_idx = curr_idx + 1;
            while(end_idx < motif_sites.size()) {
                if(motif_sites[end_idx] - motif_sites[end_idx - 1] > calling_parameters.min_separation)
                    break;
                end_idx += 1;
            }
            groups.push_back(std::make_pair(curr_idx, end_idx));
            curr_idx = end_idx;
        }

        for(size_t group_idx = 0; group_idx < groups.size(); ++group_idx) {

            size_t start_idx = groups[group_idx].first;
            size_t end_idx = groups[group_idx].second;

            // the coordinates on the reference substring for this group of sites
            int sub_start_pos = motif_sites[start_idx] - calling_parameters.min_flank;
            int sub_end_pos = motif_sites[end_idx - 1] + calling_parameters.min_flank;
            int span = motif_sites[end_idx - 1] - motif_sites[start_idx];

            // skip if too close to the start of the read alignment or
            // if the reference range is too large to efficiently call
            if(sub_start_pos <= calling_parameters.min_separation || span > 200) {
                continue;
            }

            std::string subseq = ref_seq.substr(sub_start_pos, sub_end_pos - sub_start_pos + 1);
            std::string rc_subseq = calling_parameters.alphabet->reverse_complement(subseq);

            int calling_start = sub_start_pos + ref_start_pos;
            int calling_end = sub_end_pos + ref_start_pos;

            // using the reference-to-event map, look up the event indices for this segment
            int e1,e2;
            bool bounded = AlignmentDB::_find_by_ref_bounds(event_align_record.aligned_events,
                                                            calling_start,
                                                            calling_end,
                                                            e1,
                                                            e2);

            double ratio = fabs(e2 - e1) / (calling_start - calling_end);

            // Only process this region if the the read is aligned within the boundaries
            // and the span between the start/end is not unusually short
            if(!bounded || abs(e2 - e1) <= 10 || ratio > MAX_EVENT_TO_BP_RATIO) {
                continue;
            }

            uint32_t hmm_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;

            // Set up event data
            HMMInputData data;
            data.read = &sr;
            data.pore_model = sr.get_model(strand_idx, calling_parameters.methylation_type);
            data.strand = strand_idx;
            data.rc = event_align_record.rc;
            data.event_start_idx = e1;
            data.event_stop_idx = e2;
            data.event_stride = data.event_start_idx <= data.event_stop_idx ? 1 : -1;

            // Calculate the likelihood of the unmethylated sequence
            HMMInputSequence unmethylated(subseq, rc_subseq, calling_parameters.alphabet);
            double unmethylated_score = profile_hmm_score(unmethylated, data, hmm_flags);

            // Methylate all motifs in the sequence and score again
            std::string m_subseq = calling_parameters.alphabet->methylate(subseq);
            std::string rc_m_subseq = calling_parameters.alphabet->reverse_complement(m_subseq);

            // Calculate the likelihood of the methylated sequence
            HMMInputSequence methylated(m_subseq, rc_m_subseq, calling_parameters.alphabet);
            double methylated_score = profile_hmm_score(methylated, data, hmm_flags);

            // Aggregate score
            int start_position = motif_sites[start_idx] + ref_start_pos;
            int end_position = motif_sites[end_idx - 1] + ref_start_pos;
            
            // do not output if outside the window boundaries
            if((region_start != -1 && start_position < region_start) ||
               (region_end != -1 && end_position >= region_end)) {
                continue;
            }

            auto iter = site_score_map->find(start_position);
            if(iter == site_score_map->end()) {
                // insert new score into the map
                ScoredSite ss;
                ss.chromosome = contig;
                ss.start_position = start_position;
                ss.end_position = end_position;
                ss.n_motif = end_idx - start_idx;

                // extract the motif site(s) with a k-mers worth of surrounding context
                size_t site_output_start = motif_sites[start_idx] - k + 1;
                size_t site_output_end =  motif_sites[end_idx - 1] + k;
                ss.sequence = ref_seq.substr(site_output_start, site_output_end - site_output_start);

                // insert into the map
                iter = site_score_map->insert(std::make_pair(start_position, ss)).first;
            }

            // set strand-specific score
            // upon output below the strand scores will be summed
            iter->second.ll_unmethylated[strand_idx] = unmethylated_score;
            iter->second.ll_methylated[strand_idx] = methylated_score;
            iter->second.strands_scored += 1;
        } // for group
    } // for strands

#if 0
    #pragma omp critical(call_methylation_write)
    {
        // write all sites for this read
        for(auto iter = site_score_map->begin(); iter != site_score_map->end(); ++iter) {

            const ScoredSite& ss = iter->second;
            double sum_ll_m = ss.ll_methylated[0] + ss.ll_methylated[1];
            double sum_ll_u = ss.ll_unmethylated[0] + ss.ll_unmethylated[1];
            double diff = sum_ll_m - sum_ll_u;

            // do not output if outside the window boundaries
            if((region_start != -1 && ss.start_position < region_start) ||
               (region_end != -1 && ss.end_position >= region_end)) {
                continue;
            }

            fprintf(handles.site_writer, "%s\t%s\t%d\t%d\t", ss.chromosome.c_str(), read_orientation.c_str(), ss.start_position, ss.end_position);
            fprintf(handles.site_writer, "%s\t%.2lf\t", sr.read_name.c_str(), diff);
            fprintf(handles.site_writer, "%.2lf\t%.2lf\t", sum_ll_m, sum_ll_u);
            fprintf(handles.site_writer, "%d\t%d\t%s\n", ss.strands_scored, ss.n_motif, ss.sequence.c_str());
        }
        
        // if the write bam option is turned on, write the alignment to disk
        if(handles.bam_writer != NULL) {
            bam1_t* modbam_record = create_modbam_record(record, ref_seq, ref_start_pos, *site_score_map, calling_parameters);
            
            int write_ret = sam_write1(handles.bam_writer, hdr, modbam_record);
            if(write_ret < 0) {
                fprintf(stderr, "error writing bam %d\n", write_ret);
            }

            bam_destroy1(modbam_record); // automatically frees malloc'd segment
        }

    }
#endif
}

