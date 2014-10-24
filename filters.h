#ifndef __FILTERS_H__
#define __FILTERS_H__
#include <fftw3.h>
#include <gsl/gsl_wavelet.h>
#include "common.h"

#define FFTW_NTHREADS_DEFAULT 4
#define FFTW_FLAGS_DEFAULT FFTW_ESTIMATE

typedef struct filters_handle 
{
    /* public for i/o */
    ANALYSIS_WAVEFORM_BASE_TYPE *inWav;
    ANALYSIS_WAVEFORM_BASE_TYPE *outWav;
    ANALYSIS_WAVEFORM_BASE_TYPE *respWav;
    WAVELET_BASE_TYPE *waveletWav;

    /* private for internal work */
    size_t wavLen, respLen;
    int malloced;

    /* fft */
    int fftUsed;
    size_t fftLen;
    size_t fftwNThreads;
    unsigned fftwFlags;
    FFTW(plan) fftwPlan, fftwPlan1, fftwPlan2;
    FFT_BASE_TYPE *fftwWork, *fftwWork1;
    /* window function */
    FFT_BASE_TYPE *fftwWin;
    double fftwS1, fftwS2, dt;

    /* wavelet */
    gsl_wavelet *gslDWT;
    gsl_wavelet_workspace *gslDWTWork;
} filters_t;
/* If inWav == NULL, a space is malloced and freed upon closure. 
 * If inWav is given, it is not freed.
 * n is the waveform length.  The output is also of length n */
filters_t *filters_init(ANALYSIS_WAVEFORM_BASE_TYPE *inWav, size_t n);
 /* for convolution and fft */
filters_t *filters_init_for_convolution(ANALYSIS_WAVEFORM_BASE_TYPE *inWav, size_t n, size_t np);
int filters_close(filters_t *fHdl);

int filters_SavitzkyGolay(filters_t *fHdl, int m, int ld);
int filters_raisedCosine(filters_t *fHdl);
int filters_convolute(filters_t *fHdl);
int filters_hanning_window(filters_t *fHdl);
/* compute the spectrum in Fourier space, requires init_for_convolution with np = 0 */
int filters_fft_spectrum(filters_t *fHdl);
int filters_DWT(filters_t *fHdl); /* discrete wavelet transform */

#endif // __FILTERS_H__
