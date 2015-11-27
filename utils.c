#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "common.h"
#include "utils.h"

/* Random number generator incorporated from NR3 */
static uint64_t rand_u, rand_v=4101842887655102017LL, rand_w=1LL;

void rand_init(uint64_t j)
/* j is the seed, e.g. 1237026722LL */
{
    rand_v = 4101842887655102017LL;
    rand_w = 1LL;
    
    rand_u = j ^ rand_v; rand_int64();
    rand_v = rand_u; rand_int64();
    rand_w = rand_v; rand_int64();
}

uint64_t rand_int64()
{
    uint64_t x;
    
    rand_u = rand_u * 2862933555777941757LL + 7046029254386353087LL;
    rand_v ^= rand_v >> 17; rand_v ^= rand_v << 31; rand_v ^= rand_v >> 8;
    rand_w = 4294957665U*(rand_w & 0xffffffff) + (rand_w >> 32);
    x = rand_u ^ (rand_u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + rand_v) ^ rand_w;
}

double rand0_1()
/* Uniform deviation between 0 and 1 */
{
    return 5.42101086242752217E-20 * rand_int64();
}

double rand_gauss()
/* Gauss (Deviates) Distribution */
{
    double v1, v2, R, fac;
    static double gset;
    static int iset;
    if(iset == 0) {
        do {
            v1 = 2.0 * rand0_1() - 1.0;
            v2 = 2.0 * rand0_1() - 1.0;
            R = v1*v1 + v2*v2;
        } while(R >= 1.0 || R == 0.0);
        fac = sqrt(-2.0 * log(R)/R);
        gset = v1 * fac;
        iset = 1;
        return (v2 * fac);
    } else {
        iset = 0;
        return gset;
    }
}

double rand_exp(double alpha)
/* Exponential Distribution */
{
    double u;
    do { u = rand0_1();
    } while (u == 0.0);
    return -log(u)/alpha;
}

static inline size_t qs_partition(QS_TYPE *list, size_t left, size_t right, size_t pivotIndex)
{
    QS_TYPE pivotValue, tmp;
    size_t i, storeIndex;

    pivotValue = list[pivotIndex];
    /* move pivot to the end */
    tmp = list[pivotIndex]; list[pivotIndex] = list[right]; list[right] = tmp;
    storeIndex = left;
    for(i=left; i<=right-1; i++) {
        if(list[i] < pivotValue) {
            tmp = list[storeIndex]; list[storeIndex] = list[i]; list[i] = tmp;
            storeIndex++;
        }
    }
    /* move pivot to its final place */
    tmp = list[right]; list[right] = list[storeIndex]; list[storeIndex] = tmp;
    return storeIndex;
}

static inline QS_TYPE qs_select(QS_TYPE *list, size_t left, size_t right, size_t kth)
{
    size_t pivotIndex;
    if(left == right) return list[left];
    pivotIndex = left + rand_int64() % (right - left + 1LL);
    pivotIndex = qs_partition(list, left, right, pivotIndex);
    if(kth == pivotIndex) {
        return list[kth];
    } else if (kth < pivotIndex) {
        return qs_select(list, left, pivotIndex - 1, kth);
    } else {
        return qs_select(list, pivotIndex + 1, right, kth);
    }
}

QS_TYPE quickselect(QS_TYPE *a, size_t n, size_t kth)
/* Quickselect for finding median or kth smallest element, k starts from 0 */
{
    QS_TYPE *list, v;
    list = (QS_TYPE*)malloc(sizeof(QS_TYPE) * n);
    memcpy(list, a, sizeof(QS_TYPE)*n);
    v = qs_select(list, 0, n-1, kth);
    free(list);
    return v;
}

#ifdef UTILS_DEBUG_ENABLEMAIN
int main(int argc, char **argv)
{
    size_t i, k;
    uint64_t rand_seed=1237026722LL;
    QS_TYPE list[16];

    rand_init(rand_seed);
    for(i=0; i<16; i++) {
        list[i] = (QS_TYPE)(rand_int64() % 16LL);
        printf(" %g", list[i]);
    }
    printf("\n");
    for(k = 0; k<16; k++) {
        printf("%zdth smallest element is %g\n", k, quickselect(list, 16, k));
    }
    
    return EXIT_SUCCESS;
}
#endif
