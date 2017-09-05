#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "event_detection.h"
#include "scrappie_stdlib.h"

typedef struct {
    int DEF_PEAK_POS;
    float DEF_PEAK_VAL;
    float *signal;
    size_t signal_length;
    float threshold;
    size_t window_length;
    size_t masked_to;
    int peak_pos;
    float peak_value;
    bool valid_peak;
} Detector;
typedef Detector *DetectorPtr;

/**
 *   Compute cumulative sum and sum of squares for a vector of data
 *
 *   Element i  sum (sumsq) is the sum (sum of squares) up to but
 *   excluding element i of the inputy data.
 *
 *   @param data      double[d_length]   Data to be summed over (in)
 *   @param sum       double[d_length + 1]   Vector to store sum (out)
 *   @param sumsq     double[d_length + 1]   Vector to store sum of squares (out)
 *   @param d_length                     Length of data vector
 **/
void compute_sum_sumsq(const float *data, double *sum,
                       double *sumsq, size_t d_length) {
    RETURN_NULL_IF(NULL == data, );
    RETURN_NULL_IF(NULL == sum, );
    RETURN_NULL_IF(NULL == sumsq, );
    assert(d_length > 0);

    sum[0] = 0.0f;
    sumsq[0] = 0.0f;
    for (size_t i = 0; i < d_length; ++i) {
        sum[i + 1] = sum[i] + data[i];
        sumsq[i + 1] = sumsq[i] + data[i] * data[i];
    }
}

/**
 *   Compute windowed t-statistic from summary information
 *
 *   @param sum       double[d_length]  Cumulative sums of data (in)
 *   @param sumsq     double[d_length]  Cumulative sum of squares of data (in)
 *   @param d_length                    Length of data vector
 *   @param w_length                    Window length to calculate t-statistic over
 *
 *   @returns float array containing tstats.  Returns NULL on error
 **/
float *compute_tstat(const double *sum, const double *sumsq,
                     size_t d_length, size_t w_length) {
    assert(d_length > 0);
    assert(w_length > 0);
    RETURN_NULL_IF(NULL == sum, NULL);
    RETURN_NULL_IF(NULL == sumsq, NULL);

    float *tstat = calloc(d_length, sizeof(float));
    RETURN_NULL_IF(NULL == tstat, NULL);

    const float eta = FLT_MIN;
    const float w_lengthf = (float)w_length;

    // Quick return:
    //   t-test not defined for number of points less than 2
    //   need at least as many points as twice the window length
    if (d_length < 2 * w_length || w_length < 2) {
        for (size_t i = 0; i < d_length; ++i) {
            tstat[i] = 0.0f;
        }
        return tstat;
    }
    // fudge boundaries
    for (size_t i = 0; i < w_length; ++i) {
        tstat[i] = 0;
        tstat[d_length - i - 1] = 0;
    }

    // get to work on the rest
    for (size_t i = w_length; i <= d_length - w_length; ++i) {
        double sum1 = sum[i];
        double sumsq1 = sumsq[i];
        if (i > w_length) {
            sum1 -= sum[i - w_length];
            sumsq1 -= sumsq[i - w_length];
        }
        float sum2 = (float)(sum[i + w_length] - sum[i]);
        float sumsq2 = (float)(sumsq[i + w_length] - sumsq[i]);
        float mean1 = sum1 / w_lengthf;
        float mean2 = sum2 / w_lengthf;
        float combined_var = sumsq1 / w_lengthf - mean1 * mean1
            + sumsq2 / w_lengthf - mean2 * mean2;

        // Prevent problem due to very small variances
        combined_var = fmaxf(combined_var, eta);

        //t-stat
        //  Formula is a simplified version of Student's t-statistic for the
        //  special case where there are two samples of equal size with
        //  differing variance
        const float delta_mean = mean2 - mean1;
        tstat[i] = fabs(delta_mean) / sqrt(combined_var / w_lengthf);
    }

    return tstat;
}

/**
 *
 *   @returns array of length nsample whose elements contain peak positions
 *   Remaining elements are padded by zeros.
 **/
size_t *short_long_peak_detector(DetectorPtr short_detector,
                                 DetectorPtr long_detector,
                                 const float peak_height) {
    assert(short_detector->signal_length == long_detector->signal_length);
    RETURN_NULL_IF(NULL == short_detector->signal, NULL);
    RETURN_NULL_IF(NULL == long_detector->signal, NULL);

    const size_t ndetector = 2;
    DetectorPtr detectors[] = { short_detector, long_detector };

    size_t *peaks = calloc(short_detector->signal_length, sizeof(size_t));
    RETURN_NULL_IF(NULL == peaks, NULL);

    size_t peak_count = 0;
    for (size_t i = 0; i < short_detector->signal_length; i++) {
        for (int k = 0; k < ndetector; k++) {
            DetectorPtr detector = detectors[k];
            //Carry on if we've been masked out
            if (detector->masked_to >= i) {
                continue;
            }

            float current_value = detector->signal[i];

            if (detector->peak_pos == detector->DEF_PEAK_POS) {
                //CASE 1: We've not yet recorded a maximum
                if (current_value < detector->peak_value) {
                    //Either record a deeper minimum...
                    detector->peak_value = current_value;
                } else if (current_value - detector->peak_value >
                           peak_height) {
                    // ...or we've seen a qualifying maximum
                    detector->peak_value = current_value;
                    detector->peak_pos = i;
                    //otherwise, wait to rise high enough to be considered a peak
                }
            } else {
                //CASE 2: In an existing peak, waiting to see if it is good
                if (current_value > detector->peak_value) {
                    //Update the peak
                    detector->peak_value = current_value;
                    detector->peak_pos = i;
                }
                //Dominate other tstat signals if we're going to fire at some point
                if (detector == short_detector) {
                    if (detector->peak_value > detector->threshold) {
                        long_detector->masked_to =
                            detector->peak_pos + detector->window_length;
                        long_detector->peak_pos =
                            long_detector->DEF_PEAK_POS;
                        long_detector->peak_value =
                            long_detector->DEF_PEAK_VAL;
                        long_detector->valid_peak = false;
                    }
                }
                //Have we convinced ourselves we've seen a peak
                if (detector->peak_value - current_value > peak_height
                    && detector->peak_value > detector->threshold) {
                    detector->valid_peak = true;
                }
                //Finally, check the distance if this is a good peak
                if (detector->valid_peak
                    && (i - detector->peak_pos) >
                    detector->window_length / 2) {
                    //Emit the boundary and reset
                    peaks[peak_count] = detector->peak_pos;
                    peak_count++;
                    detector->peak_pos = detector->DEF_PEAK_POS;
                    detector->peak_value = current_value;
                    detector->valid_peak = false;
                }
            }
        }
    }

    return peaks;
}

/**  Create an event given boundaries
 *
 *   Note: Bounds are CADLAG (i.e. lower bound is contained in the interval but
 *   the upper bound is not).
 *
 *  @param start Index of lower bound
 *  @param end Index of upper bound
 *  @param sums
 *  @param sumsqs
 *  @param nsample  Total number of samples in read
 *
 *  @returns An initialised event.  A 'null' event is returned on error.
 **/
event_t create_event(size_t start, size_t end, double const *sums,
                     double const *sumsqs, size_t nsample) {
    assert(start < nsample);
    assert(end <= nsample);

    event_t event = { 0 };
    event.pos = -1;
    event.state = -1;
    RETURN_NULL_IF(NULL == sums, event);
    RETURN_NULL_IF(NULL == sumsqs, event);

    event.start = (uint64_t)start;
    event.length = (float)(end - start);
    event.mean = (float)(sums[end] - sums[start]) / event.length;
    const float deltasqr = (sumsqs[end] - sumsqs[start]);
    const float var = deltasqr / event.length - event.mean * event.mean;
    event.stdv = sqrtf(fmaxf(var, 0.0f));

    return event;
}

event_table create_events(size_t const *peaks, double const *sums,
                          double const *sumsqs, size_t nsample) {
    event_table et = { 0 };
    RETURN_NULL_IF(NULL == sums, et);
    RETURN_NULL_IF(NULL == sumsqs, et);
    RETURN_NULL_IF(NULL == peaks, et);

    // Count number of events found
    size_t n = 1;
    for (size_t i = 0; i < nsample; ++i) {
        if (peaks[i] > 0 && peaks[i] < nsample) {
            n++;
        }
    }

    et.event = calloc(n, sizeof(event_t));
    RETURN_NULL_IF(NULL == et.event, et);

    et.n = n;
    et.end = et.n;


    // First event -- starts at zero
    et.event[0] = create_event(0, peaks[0], sums, sumsqs, nsample);
    // Other events -- peak[i-1] -> peak[i]
    for(size_t ev=1 ; ev < n - 1 ; ev++){
        et.event[ev] = create_event(peaks[ev - 1], peaks[ev], sums, sumsqs, nsample);
    }
    // Last event -- ends at nsample
    et.event[n - 1] = create_event(peaks[n - 2], nsample, sums, sumsqs, nsample);

    return et;
}

event_table detect_events(raw_table const rt, detector_param const edparam) {

    event_table et = { 0 };
    RETURN_NULL_IF(NULL == rt.raw, et);

    double *sums = calloc(rt.n + 1, sizeof(double));
    double *sumsqs = calloc(rt.n + 1, sizeof(double));

    compute_sum_sumsq(rt.raw, sums, sumsqs, rt.n);
    float *tstat1 = compute_tstat(sums, sumsqs, rt.n, edparam.window_length1);
    float *tstat2 = compute_tstat(sums, sumsqs, rt.n, edparam.window_length2);

    Detector short_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = FLT_MAX,
        .signal = tstat1,
        .signal_length = rt.n,
        .threshold = edparam.threshold1,
        .window_length = edparam.window_length1,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = FLT_MAX,
        .valid_peak = false
    };

    Detector long_detector = {
        .DEF_PEAK_POS = -1,
        .DEF_PEAK_VAL = FLT_MAX,
        .signal = tstat2,
        .signal_length = rt.n,
        .threshold = edparam.threshold2,
        .window_length = edparam.window_length2,
        .masked_to = 0,
        .peak_pos = -1,
        .peak_value = FLT_MAX,
        .valid_peak = false
    };

    size_t *peaks =
        short_long_peak_detector(&short_detector, &long_detector,
                                 edparam.peak_height);

    et = create_events(peaks, sums, sumsqs, rt.n);

    free(peaks);
    free(tstat2);
    free(tstat1);
    free(sumsqs);
    free(sums);

    return et;
}
