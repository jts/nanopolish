//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_scorereads -- score reads against an alignment, model
//

#ifndef NANOPOLISH_SCOREREADS_H
#define NANOPOLISH_SCOREREADS_H

double model_score(SquiggleRead &sr,
                   const size_t strand_idx,
                   const faidx_t *fai, 
                   const std::vector<EventAlignment> &alignment_output,
                   const size_t events_per_segment);

int scorereads_main(int argc, char** argv);

#endif 
