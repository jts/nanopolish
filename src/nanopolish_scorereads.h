//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_scorereads -- score reads against an alignment, model
//

#ifndef NANOPOLISH_SCOREREADS_H
#define NANOPOLISH_SCOREREADS_H


std::vector<EventAlignment> alignment_from_read(SquiggleRead& sr,
                                                const size_t strand_idx,
                                                const size_t read_idx,
                                                const std::string& alternative_model_type,
                                                const faidx_t* fai,
                                                const bam_hdr_t* hdr,
                                                const bam1_t* record,
                                                int region_start,
                                                int region_end);

double model_score(SquiggleRead &sr,
                   const size_t strand_idx,
                   const faidx_t *fai, 
                   const std::vector<EventAlignment> &alignment_output,
                   const size_t events_per_segment,
                   TransitionParameters* transition_training);

int scorereads_main(int argc, char** argv);

#endif 
