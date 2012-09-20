#ifndef __FILTERS_H__
#define __FILTERS_H__
#include <fftw3.h>
#include <gsl/gsl_wavelet.h>
#include "common.h"

typedef struct filters_handle 
{
    /* public for i/o */
    ANALYSIS_WAVEFORM_BASE_TYPE *inWav;
    ANALYSIS_WAVEFORM_BASE_TYPE *outWav;
    ANALYSIS_WAVEFORM_BASE_TYPE *respWav;
    WAVELET_BASE_TYPE *waveletWav;

    /* private for internal work */
    size_t wavLen, respLen;
    int fftUsed;
    size_t fftLen;
    FFTW(plan) fftwPlan, fftwPlan1, fftwPlan2;
    FFT_BASE_TYPE *fftWork, *fftWork1;

    gsl_wavelet *gslDWT;
    gsl_wavelet_workspace *gslDWTWork;
} filters_t;

filters_t *filters_init(size_t n); /* input the waveform length.  The output is also of length n */
filters_t *filters_init_for_convolution(size_t n, size_t np); /* for convolution */
int filters_close(filters_t *fHdl);

int filters_SavitzkyGolay(filters_t *fHdl, int m, int ld);
int filters_raisedCosine(filters_t *fHdl);
int filters_convolute(filters_t *fHdl);
int filters_DWT(filters_t *fHdl); /* discrete wavelet transform */

#endif // __FILTERS_H__
