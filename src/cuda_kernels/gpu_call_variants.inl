#include <cuda_kernels/GpuAligner.h>
#include <thread>
#include <chrono>
#include <future>

void prepareForBaseEditCandidates(int start,
                                  int end,
                                  const AlignmentDB& alignments,
                                  std::string contig,
                                  std::vector<std::vector<Variant>> &tmp_variants_vector,
                                  std::vector<Haplotype> &haplotypes,
                                  std::vector<std::vector<HMMInputData>> &event_sequences_vector){
    for(int i = start; i<=end; i++){
        int calling_start = i - opt::screen_flanking_sequence;
        int calling_end = i + 1 + opt::screen_flanking_sequence;

        if (!alignments.are_coordinates_valid(contig, calling_start, calling_end)) {
            return;
        }

        std::vector<Variant> tmp_variants;
        for (size_t j = 0; j < 4; ++j) {
            // Substitutions
            Variant v;
            v.ref_name = contig;
            v.ref_position = i;
            v.ref_seq = alignments.get_reference_substring(contig, i, i);
            v.alt_seq = "ACGT"[j];

            if (v.ref_seq != v.alt_seq) {
                tmp_variants.push_back(v);
            }

            // Insertions
            v.alt_seq = v.ref_seq + "ACGT"[j];
            // ignore insertions of the type "A" -> "AA" as these are redundant
            if (v.alt_seq[1] != v.ref_seq[0]) {
                tmp_variants.push_back(v);
            }
        }

        // deletion
        Variant del;
        del.ref_name = contig;
        del.ref_position = i - 1;
        del.ref_seq = alignments.get_reference_substring(contig, i - 1, i);
        del.alt_seq = del.ref_seq[0];

        // ignore deletions of the type "AA" -> "A" as these are redundant
        if (del.alt_seq[0] != del.ref_seq[1]) {
            tmp_variants.push_back(del);
        }

        // Screen variants by score
        // We do this internally here as it is much faster to get the event sequences
        // for the entire window for all variants at this position once, rather than
        // for each variant individually
        std::vector<HMMInputData> event_sequences = alignments.get_event_subsequences(contig, calling_start, calling_end);

        Haplotype test_haplotype(contig,
                                 calling_start,
                                 alignments.get_reference_substring(contig,
                                                                    calling_start,
                                                                    calling_end));

        haplotypes.push_back(test_haplotype);
        event_sequences_vector.push_back(event_sequences);
        tmp_variants_vector.push_back(tmp_variants);
    }
}


void locusRangeBaseEditCandidateGPU(int start,
                                    int end,
                                    const AlignmentDB& alignments,
                                    uint32_t alignment_flags,
                                    std::vector<Variant> &out_variants,
                                    std::string contig,
                                    GpuAligner &aligner,
                                    std::mutex &outVariantsMutex) {
    std::vector<std::vector<Variant>> tmp_variants_vector;
    std::vector<Haplotype> haplotypes;
    std::vector<std::vector<HMMInputData>> event_sequences_vector;

    prepareForBaseEditCandidates(start,
                                 end,
                                 alignments,
                                 contig,
                                 tmp_variants_vector,
                                 haplotypes,
                                 event_sequences_vector);

    std::vector<Variant> scoredVariants = aligner.variantScoresThresholded(tmp_variants_vector,
                                                                           haplotypes,
                                                                           event_sequences_vector,
                                                                           alignment_flags,
                                                                           opt::screen_score_threshold,
                                                                           opt::methylation_types);
    for (auto variant: scoredVariants) {
        if (variant.quality > 0) {
            std::lock_guard<std::mutex> lock(outVariantsMutex);
            out_variants.push_back(variant);
        }
    }

}

std::vector<Variant> generate_candidate_single_base_edits_gpu(const AlignmentDB& alignments,
                                                              int region_start,
                                                              int region_end,
                                                              uint32_t alignment_flags){

    std::mutex outVariantsMutex;
    std::vector<Variant> out_variants;
    std::string contig = alignments.get_region_contig();

    // Add all positively-scoring single-base changes into the candidate set
    size_t num_workers = (opt::num_threads < MAX_NUM_WORKERS) ? opt::num_threads : MAX_NUM_WORKERS;
    std::vector<GpuAligner> gpuAligners(num_workers);

    //std::vector<std::thread> workerThreads(num_workers);
    std::vector<std::future<void>> handles(num_workers);

    int nextLocusBegin = region_start;
    int nextLocusEnd = nextLocusBegin + LOCI_PER_WORKER;
    bool finished = false;

    //Initialise the workers
    for (int workerIdx = 0; workerIdx < num_workers; workerIdx++) {
        auto aligner = std::ref(gpuAligners[workerIdx]);
        if (!finished) {
            if (nextLocusEnd == region_end) {
                finished = true;
            }
            handles[workerIdx] = std::async(std::launch::async,
                                            locusRangeBaseEditCandidateGPU,
                                            nextLocusBegin,
                                            nextLocusEnd,
                                            std::ref(alignments),
                                            alignment_flags,
                                            std::ref(out_variants),
                                            std::ref(contig),
                                            aligner,
                                            std::ref(outVariantsMutex));
            if ((nextLocusEnd + LOCI_PER_WORKER) < region_end){
                nextLocusBegin = nextLocusEnd + 1;
                nextLocusEnd = nextLocusBegin + LOCI_PER_WORKER - 1;
            }else{
                nextLocusBegin = nextLocusEnd + 1;
                nextLocusEnd = region_end;
            }
        }
    }

    //Round robin - assigning work to the workers until out of candidates
    while (!finished) {
        for (int i = 0; i < num_workers; i++) {
            auto status = handles[i].wait_for(std::chrono::microseconds(100));
            if (status == std::future_status::ready && (!finished)) {
                if (nextLocusEnd == region_end){
                    finished = true;
                }
                auto aligner = std::ref(gpuAligners[i]);
                handles[i].get();
                handles[i] = std::async(std::launch::async,
                                        locusRangeBaseEditCandidateGPU,
                                        nextLocusBegin,
                                        nextLocusEnd,
                                        std::ref(alignments),
                                        alignment_flags,
                                        std::ref(out_variants),
                                        std::ref(contig),
                                        aligner,
                                        std::ref(outVariantsMutex));
                if ((nextLocusEnd + LOCI_PER_WORKER) < region_end){
                    nextLocusBegin = nextLocusEnd + 1;
                    nextLocusEnd = nextLocusBegin + LOCI_PER_WORKER - 1;
                }else{
                    nextLocusBegin = nextLocusEnd + 1;
                    nextLocusEnd = region_end;
                }
            }
        }
    }

    //Block until all workers are complete
    for (int workerIdx = 0; workerIdx < num_workers; workerIdx++) {
        handles[workerIdx].wait();
    }
    return  out_variants;
}
