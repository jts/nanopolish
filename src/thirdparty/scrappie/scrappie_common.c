#include "scrappie_common.h"
#include "scrappie_stdlib.h"
//#include "util.h"
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

int floatcmp(const void *x, const void *y) {
    float d = *(float *)x - *(float *)y;
    if (d > 0) {
        return 1;
    }
    return -1;
}

/**  Quantiles from n array
 *
 *  Using a relatively inefficent qsort resulting in O(n log n)
 *  performance but better performance is possible for small np.
 *  The array p is modified inplace, containing which quantiles to
 *  calculation on input and the quantiles on output; on error, p
 *  is filled with the value NAN.
 *
 *  @param x An array to calculate quantiles from
 *  @param nx Length of array x
 *  @param p An array containing quantiles to calculate [in/out]
 *  @param np Length of array p
 *
 *  @return void
 **/
void quantilef(const float *x, size_t nx, float *p, size_t np) {
    if (NULL == p) {
        return;
    }
    for (int i = 0; i < np; i++) {
        assert(p[i] >= 0.0f && p[i] <= 1.0f);
    }
    if (NULL == x) {
        for (int i = 0; i < np; i++) {
            p[i] = NAN;
        }
        return;
    }
    // Sort array
    float *space = malloc(nx * sizeof(float));
    if (NULL == space) {
        for (int i = 0; i < np; i++) {
            p[i] = NAN;
        }
        return;
    }
    memcpy(space, x, nx * sizeof(float));
    qsort(space, nx, sizeof(float), floatcmp);

    // Extract quantiles
    for (int i = 0; i < np; i++) {
        const size_t idx = p[i] * (nx - 1);
        const float remf = p[i] * (nx - 1) - idx;
        if (idx < nx - 1) {
            p[i] = (1.0 - remf) * space[idx] + remf * space[idx + 1];
        } else {
            // Should only occur when p is exactly 1.0
            p[i] = space[idx];
        }
    }

    free(space);
    return;
}

/** Median of an array
 *
 *  Using a relatively inefficent qsort resulting in O(n log n)
 *  performance but O(n) is possible.
 *
 *  @param x An array to calculate median of
 *  @param n Length of array
 *
 *  @return Median of array on success, NAN otherwise.
 **/
float medianf(const float *x, size_t n) {
    float p = 0.5;
    quantilef(x, n, &p, 1);
    return p;
}

/** Median Absolute Deviation of an array
 *
 *  @param x An array to calculate the MAD of
 *  @param n Length of array
 *  @param med Median of the array.  If NAN then median is calculated.
 *
 *  @return MAD of array on success, NAN otherwise.
 **/
float madf(const float *x, size_t n, const float *med) {
    const float mad_scaling_factor = 1.4826;
    if (NULL == x) {
        return NAN;
    }
    if (1 == n) {
        return 0.0f;
    }

    float *absdiff = malloc(n * sizeof(float));
    if (NULL == absdiff) {
        return NAN;
    }

    const float _med = (NULL == med) ? medianf(x, n) : *med;

    for (size_t i = 0; i < n; i++) {
        absdiff[i] = fabsf(x[i] - _med);
    }

    const float mad = medianf(absdiff, n);
    free(absdiff);
    return mad * mad_scaling_factor;
}


raw_table trim_and_segment_raw(raw_table rt, int trim_start, int trim_end, int varseg_chunk, float varseg_thresh) {
    RETURN_NULL_IF(NULL == rt.raw, (raw_table){0});

    rt = trim_raw_by_mad(rt, varseg_chunk, varseg_thresh);
    RETURN_NULL_IF(NULL == rt.raw, (raw_table){0});

    rt.start += trim_start;
    rt.end -= trim_end;

    if (rt.start >= rt.end) {
        free(rt.raw);
        return (raw_table){0};
    }

    return rt;
}

/**  Simple segmentation of a raw read by thresholding the MAD
 *
 *  The MAD of the raw signal is calculated for non-overlapping chunks and then
 *  thresholded to find regions at the beginning and end of the signal that have
 *  unusually low variation (generally a stall or open pore).  The threshhold is
 *  derived from the distribution of the calaculated MADs.
 *
 *  The threshold is chosen to be high since a single chunk above it will trigger
 *  the end of the trimming: the threshhold is chosen so it is unlikely to be
 *  exceeded in the leader but commonly exceeded in the main read.
 *
 *  @param rt Structure containing raw signal
 *  @param chunk_size Size of non-overlapping chunks
 *  @param perc  The quantile to be calculated to use for threshholding
 *
 *  @return A range structure containing new start and end for read
 **/
raw_table trim_raw_by_mad(raw_table rt, int chunk_size, float perc) {
    assert(chunk_size > 1);
    assert(perc >= 0.0 && perc <= 1.0);

    const size_t nsample = rt.end - rt.start;
    const size_t nchunk = nsample / chunk_size;
    // Truncation of end to be consistent with Sloika
    rt.end = nchunk * chunk_size;

    float *madarr = malloc(nchunk * sizeof(float));
    RETURN_NULL_IF(NULL == madarr, (raw_table){0});
    for (size_t i = 0; i < nchunk; i++) {
        madarr[i] = madf(rt.raw + rt.start + i * chunk_size, chunk_size, NULL);
    }
    quantilef(madarr, nchunk, &perc, 1);

    const float thresh = perc;
    for (size_t i = 0; i < nchunk; i++) {
        if (madarr[i] > thresh) {
            break;
        }
        rt.start += chunk_size;
    }
    for (size_t i = nchunk; i > 0; i--) {
        if (madarr[i - 1] > thresh) {
            break;
        }
        rt.end -= chunk_size;
    }
    assert(rt.end > rt.start);

    free(madarr);

    return rt;
}
