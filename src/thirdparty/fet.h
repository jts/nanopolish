#ifndef FET_H
#define FET_H

#ifdef __cplusplus
extern "C" {
#endif

double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);

#ifdef __cplusplus
}
#endif

#endif
