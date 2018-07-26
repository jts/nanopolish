#include <vector>
#include "nanopolish_variant.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <algorithm>
#include <queue>
#include <sstream>
#include <fstream>
#include <set>
#include <omp.h>
#include <getopt.h>
#include <iterator>
#include "htslib/faidx.h"
#include "nanopolish_poremodel.h"
#include "nanopolish_transition_parameters.h"
#include "nanopolish_matrix.h"
#include "nanopolish_klcs.h"
#include "nanopolish_profile_hmm.h"
#include "nanopolish_alignment_db.h"
#include "nanopolish_anchor.h"
#include "nanopolish_variant.h"
#include "nanopolish_haplotype.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_duration_model.h"
#include "nanopolish_variant_db.h"
#include "profiler.h"
#include "progress.h"
#include "stdaln.h"
#include <chrono>
#include <cuda.h>
#include <cuda_runtime.h>

#ifndef GPU_ALIGNER_H
#define GPU_ALIGNER_H

#define LOCI_PER_WORKER 32
#define MAX_COVERAGE 300
#define MAX_SEQUENCE_LENGTH 50
#define MAX_NUM_VARIANTS_PER_LOCUS 10

//Data to be scored
typedef struct {
    std::vector<HMMInputSequence> stateSequences;
    std::vector<HMMInputData> rawData;
} ScoreSet;

class GpuAligner
{
public:
    GpuAligner();
    ~GpuAligner();

    std::vector<Variant>
      variantScoresThresholded(std::vector<std::vector<Variant>>,
			       std::vector<Haplotype>,
			       std::vector<std::vector<HMMInputData>>,
              uint32_t alignment_flags, int screen_score_threshold, std::vector<std::string> methylation_types);

    std::vector<std::vector<double>> scoreKernel(std::vector<HMMInputSequence> sequences,
    std::vector<HMMInputData> event_sequences,
            uint32_t alignment_flags);
    std::vector<std::vector<std::vector<double>>> scoreKernelMod(std::vector<ScoreSet> &scoreSets,
                                                                             uint32_t alignment_flags);
private:
    float* scaleDev;
    float* shiftDev;
    float* varDev;
    float* logVarDev;
    float * eventMeans;
    float * preFlankingHost;
    float * postFlankingHost;
    int* eventOffsetsDev;
    int* eventStridesDev;
    int* eventStartsDev;
    int* numRowsDev;
    float* postFlankingDev;
    float* preFlankingDev;
    float* eventMeansDev;
    float* eventsPerBaseDev;
    float* poreModelLevelStdvDev;
    float* poreModelLevelLogStdvDev;
    float* poreModelLevelMeanDev;
    float* poreModelDev;
    float* poreModelHost;
    int * sequenceOffsetsDev;

    // NEW - for MOD kernel
    int * readLengthsHost;
    int * eventStartsHost;
    int * eventStridesHost;
    float * eventsPerBaseHost;
    float * scaleHost;
    float * shiftHost;
    float * varHost;
    float * logVarHost;
    int * sequenceLengthsHost;
    int * eventOffsetsHost;
    int * sequenceOffsetsHost;
    int * readIdxHost;
    int * seqIdxHost;

    int * readLengthsDev;
    int * sequenceLengthsDev;
    int * readIdxDev;
    int * seqIdxDev;

    float * returnValuesHost;
    float * scoresDev;

    int * kmerRanks;
    int * kmerRanksDev;

    bool poreModelInitialized;

    std::vector<int*> kmerRanksDevPointers;
    std::vector<float*> returnValuesDevResultsPointers;
    std::vector<float*> returnValuesHostResultsPointers;

    cudaStream_t streams[8]; // TODO 8 should not be hardcoded here
};

#endif // GPU_ALIGNER_H
