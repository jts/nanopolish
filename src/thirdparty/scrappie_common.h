#pragma once
#ifndef SCRAPPIE_COMMON_H
#define SCRAPPIE_COMMON_H

#include "scrappie_structures.h"

raw_table trim_and_segment_raw(raw_table rt, int trim_start, int trim_end, int varseg_chunk, float varseg_thresh);
raw_table trim_raw_by_mad(raw_table rt, int chunk_size, float proportion);

#endif /* SCRAPPIE_COMMON_H */
