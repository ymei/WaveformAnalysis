/** \file utils.h
 * Utilities such as random number/distribution generators.
 */
#ifndef __UTILS_H__
#define __UTILS_H__
#include <stdint.h>

#define QS_TYPE double

/** Random number generator incorporated from NR3
 * @param[in] j the seed, e.g. 1237026722LL
 */
void rand_init(uint64_t j);
uint64_t rand_int64();
double rand0_1();
double rand_gauss();
/** n-dimmensional Gaussian (multivariate normal) distribution
 * @param[in] n dimension.
 * @param[in] L Cholesky decomposed covariance matrix.
 * @param[in] mean vector of mean values.
 * @param[out] pt n-dimmensional Gaussian deviate output.  pt has to be pre-allocated.
 */
double *rand_gaussnd(size_t n, const double *L, const double *mean, double *pt);
/** Exponential-decay distribution 
 * @param[in] alpha \f$\exp(-\alpha t)\f$
 */
double rand_exp(double alpha);
/** Quickselect for finding median or kth smallest element, k starts from 0.
 * @param[in] a input array, will be destroyed.
 * @param[in] n length of array a.
 * @param[in] kth kth smallest element.
 * @return value of the kth smallest element in array a.
 */
QS_TYPE quickselect(QS_TYPE *a, size_t n, size_t kth);

/** Cholesky decomposition.
 * @param[in] a input array, must be square, symmetric and positive-definite, with n*n elements.
 * @param[in] n dimension of array a and L.
 * @param[out] L decomposition result, a lower triangular matrix.
 *               If (*L)=NULL is supplied, *L is allocated.
 * @return 1 when success, 0 if a is not positive-definite.
 */
int cholesky_decomp(const double *a, size_t n, double **L);

#endif /* __UTILS_H__ */
