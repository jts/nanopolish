#pragma once
#ifndef UTIL_H
#    define UTIL_H
#    include <immintrin.h>
#    include <math.h>
#    include <stdbool.h>
#    include <stdint.h>
#    include <stdio.h>
#    include "sse_mathfun.h"

#    ifdef FAST_LOG
#        define LOGFV fast_logfv
#    else
#        define LOGFV logfv
#    endif

#    ifdef FAST_EXP
#        define EXPFV fast_expfv
#    else
#        define EXPFV expfv
#    endif

#    ifdef FAST_TANH
#        define TANHFV fast_tanhfv
#    else
#        define TANHFV tanhfv
#    endif

#    ifdef FAST_LOGISTIC
#        define LOGISTICFV fast_logisticfv
#    else
#        define LOGISTICFV logisticfv
#    endif

#    ifdef FAST_ELU
#        define ELUFV fast_elufv
#    else
#        define ELUFV elufv
#    endif

/* Create a vector of  ones.  */
extern __inline __m128 __attribute__ ((__gnu_inline__, __always_inline__))
    _mm_setone_ps(void) {
    return __extension__(__m128) {
    1.0f, 1.0f, 1.0f, 1.0f};
}

extern __inline __m128d __attribute__ ((__gnu_inline__, __always_inline__))
    _mm_setone_pd(void) {
    return __extension__(__m128d) {
    1.0, 1.0};
}


/**
 *   Standard implementations
 **/
static inline float logisticf(float x) {
    return 1.0f / (1.0f + expf(-x));
}

static inline float eluf(float x){
    return (x >= 0.0f) ? x : expm1f(x);
}

// Constants for fast exp approximation.  See Schraudolph (1999)
#    define _A 12102203.161561485f
//  Minimum maximum relative error
//#    define _B 1064986822.5027076f
//  No bias at zero
#    define _B 1065353216.0f
//  Minimum RMS relative error
//#    define _B 1064866803.6193156f
//  Minimum mean relative error
//#    define _B 1064807268.0425934f
#    define _BOUND 88.02969193111305
static inline float fast_expf(float x) {
    x = fmaxf(-_BOUND, fminf(_BOUND, x));
    union {
        uint32_t i;
        float f;
    } value = {
    .i = (uint32_t) (_A * x + _B)};
    return value.f;
}

static inline float fast_logisticf(float x) {
    return 1.0 / (1.0 + fast_expf(-x));
}

static inline float fast_tanhf(float x) {
    const float y = fast_logisticf(x + x);
    return y + y - 1.0;
}

static inline float fast_eluf(float x){
    return (x >= 0.0f) ? x : (fast_expf(x) - 1.0);
}

/**
 *   Vectorised functions
 **/
static inline __m128 __attribute__ ((__always_inline__)) expfv(__m128 x) {
    __v4sf y = (__v4sf) x;
    return (__m128) exp_ps(y);
}

static inline __m128 __attribute__ ((__always_inline__)) logfv(__m128 x) {
    __v4sf y = (__v4sf) x;
    return (__m128) log_ps(y);
}

static inline __m128 __attribute__ ((__always_inline__)) logisticfv(__m128 x) {
    const __m128 ones = _mm_setone_ps();
    return _mm_div_ps(ones, _mm_add_ps(ones, expfv(-x)));
}

static inline __m128 __attribute__ ((__always_inline__)) tanhfv(__m128 x) {
    const __m128 y = logisticfv(_mm_add_ps(x, x));
    return _mm_sub_ps(_mm_add_ps(y, y), _mm_setone_ps());
}

static inline __m128 __attribute__ ((__always_inline__)) elufv(__m128 x) {
    if(0 == _mm_movemask_ps(x)){
        // All positive, early return.
        return x;
    }
    const __m128 mask = _mm_cmpge_ps(x, _mm_setzero_ps());
    const __m128 y = expfv(x) - _mm_setone_ps();
    return _mm_or_ps(_mm_and_ps(mask, x),  _mm_andnot_ps(mask, y));
}

/**
 *    Fast vectorised approximations
 **/
static inline __m128 fast_expfv(__m128 x) {
    const __m128 a = (__m128) (__v4sf) { _A, _A, _A, _A };
    const __m128 b = (__m128) (__v4sf) { _B, _B, _B, _B };
    const __m128 _bound = (__m128) (__v4sf) { _BOUND, _BOUND, _BOUND, _BOUND };
    x = _mm_max_ps(-_bound, _mm_min_ps(_bound, x));

    __m128 y = a * x + b;
    return _mm_castsi128_ps(_mm_cvtps_epi32(y));
}

static inline __m128 fast_logfv(__m128 x) {
#    define _Alogfv 8.262958294867817e-08f
#    define _Blogfv 1064872507.1541044f
    const __m128 a = (__m128) (__v4sf) { _Alogfv, _Alogfv, _Alogfv, _Alogfv };
    const __m128 b = (__m128) (__v4sf) { _Blogfv, _Blogfv, _Blogfv, _Blogfv };
    x = _mm_cvtepi32_ps(_mm_castps_si128(x));
    return a * (x - b);
}

static inline __m128 __attribute__ ((__always_inline__)) fast_logisticfv(__m128 x) {
    return _mm_rcp_ps(_mm_add_ps(_mm_setone_ps(), fast_expfv(-x)));
}

static inline __m128 __attribute__ ((__always_inline__)) fast_tanhfv(__m128 x) {
    const __m128 y = fast_logisticfv(_mm_add_ps(x, x));
    return _mm_sub_ps(_mm_add_ps(y, y), _mm_setone_ps());
}

static inline __m128 __attribute__ ((__always_inline__)) fast_elufv(__m128 x) {
    if(0 == _mm_movemask_ps(x)){
        // All positive, early return.
        return x;
    }
    const __m128 mask = _mm_cmpge_ps(x, _mm_setzero_ps());
    const __m128 y = fast_expfv(x) - _mm_setone_ps();
    return _mm_or_ps(_mm_and_ps(mask, x),  _mm_andnot_ps(mask, y));
}

int argmaxf(const float *x, int n);
int argminf(const float *x, int n);
float valmaxf(const float *x, int n);
float valminf(const float *x, int n);

static inline int iceil(int x, int y) {
    return (x + y - 1) / y;
}

static inline int ifloor(int x, int y) {
    return x / y;
}

void quantilef(const float *x, size_t nx, float *p, size_t np);
float medianf(const float *x, size_t n);
float madf(const float *x, size_t n, const float *med);
void medmad_normalise_array(float *x, size_t n);
void studentise_array_kahan(float *x, size_t n);

bool equality_array(double const * x, double const * y, size_t n, double const tol);
bool equality_arrayf(float const * x, float const * y, size_t n, float const tol);
bool equality_arrayi(int const * x, int const * y, size_t n);

#endif                          /* UTIL_H */
