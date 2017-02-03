//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_bam_processor -- framework for iterating
// over a bam file and running an arbitrary function
// on each aligned read in parallel
//
#include <assert.h>
#include <string.h>
#include "nanopolish_bam_utils.h"

void write_bam_vardata(bam1_t* record,
                      const std::string& qname,
                      const std::vector<uint32_t> cigar,
                      const std::string& seq,
                      const std::string& qual,
                      size_t aux_reserve)
{
    record->core.l_qname = qname.length() + 1; // must be null-terminated
    record->core.l_qseq = seq.size();
    record->core.n_cigar = cigar.size();
    size_t seq_len_bytes = (record->core.l_qseq + 1) / 2; // 4-bits per symbol
    
    size_t qual_length = qual.size();
    assert( (qual_length == 1 && qual == "*") || qual_length == record->core.l_qseq);
    
    // calculate length of incoming data, in bytes
    record->m_data = record->core.l_qname + // query name
                     record->core.n_cigar * 4 + // 4 bytes per cigar op
                     seq_len_bytes + // 4 bits per base for the query seq
                     qual_length + // 1 byte per quality character
                     aux_reserve;

    // nothing copied yet
    record->l_data = 0;
    
    // allocate data
    record->data = (uint8_t*)malloc(record->m_data);

    // copy q name
    assert(record->core.l_qname <= record->m_data);
    strncpy(bam_get_qname(record), 
            qname.c_str(),
            record->core.l_qname);
    record->l_data += record->core.l_qname;
    
    // copy cigar
    assert(record->l_data + record->core.n_cigar * 4 <= record->m_data);
    memcpy(bam_get_cigar(record), 
           &cigar[0],
           record->core.n_cigar * 4);
    record->l_data += record->core.n_cigar * 4;

    // copy seq
    assert(record->l_data + seq.size() / 2 <= record->m_data);
    uint8_t* seq_data = bam_get_seq(record);
    for(size_t i = 0; i < seq_len_bytes; ++i) {
        seq_data[i] = 0;
    }
    for(size_t i = 0; i < seq.size(); ++i) {
        assert(i / 2 < seq_len_bytes);
        seq_data[i / 2] |= seq_nt16_table[(int)seq[i]] << ((~i & 1) << 2);
        record->l_data += (~i & 1);
    }

    // copy quale
    assert(record->l_data + qual_length <= record->m_data);
    uint8_t* qual_data = bam_get_qual(record);
    if(qual == "*") {
        qual_data[0] = 0xff;
        record->l_data += 1;
    } else {
        for(size_t i = 0; i < qual.size(); ++i) {
            qual_data[i] = (int8_t)qual[i];
            record->l_data += 1;
        }
    }
    assert(record->l_data <= record->m_data);
}
