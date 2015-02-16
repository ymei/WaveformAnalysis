#ifndef __UTILS_H__
#define __UTILS_H__
#include <stdint.h>

#define QS_TYPE double

/* Random number generator incorporated from NR3 */
void rand_init(uint64_t j);
uint64_t rand_int64();
double rand0_1();
double rand_gauss();
double rand_exp(double alpha);
/* Quickselect for finding median or kth smallest element, k starts from 0 */
QS_TYPE quickselect(QS_TYPE *a, size_t n, size_t kth);

#endif /* __UTILS_H__ */
