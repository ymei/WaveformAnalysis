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

#endif /* __UTILS_H__ */
