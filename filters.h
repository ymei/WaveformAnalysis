/** \file filters.h
 * A collection of filters and FFT related routines.
 */
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

/** Initialize the filter.
 * @param[in] inWav waveform input.
 * If inWav == NULL, a space is malloced and freed upon closure. 
 * If inWav is given, it is not freed.
 * @param[in] n waveform length.
 * The output is also of length n.  n can be odd in principle,
 * but an even number is preferred.
 */
filters_t *filters_init(ANALYSIS_WAVEFORM_BASE_TYPE *inWav, size_t n);
/** Initialize for convolution and fft
 * @param[in] np length of response function.  Must be an odd number.
 * np=0 is reserved for fft spectrum calculations.
 */
filters_t *filters_init_for_convolution(ANALYSIS_WAVEFORM_BASE_TYPE *inWav, size_t n, size_t np);
/** Destroy the filter */
int filters_close(filters_t *fHdl);
/** Savitzky-Golay filter.
 * @param[in] m order of polynomial.
 * @param[in] np number of points (response length)
 * @param[in] ld degree of derivative
 */
int filters_SavitzkyGolay(filters_t *fHdl, int m, int ld);
/** Raised cosine as convolution kernel.
 * Outermost bin, i=(fHdl->respLen-1)/2, is set to 0.0
 **/
int filters_raisedCosine(filters_t *fHdl);
int filters_convolute(filters_t *fHdl);
int filters_hanning_window(filters_t *fHdl);
/** Compute the spectrum in Fourier space.  Requires init_for_convolution with np = 0. 
 * Computed power spectrum is stored in fftwWork(s), normalized.
 * This function handles both odd and even number of input points.
 * The resulting power spectrum has the length (int)(n/2)+1 for both even and odd n.
 * Spectra density (linearized) is stored in fftwWork,
 * Spectrum (linearized) is stored in fftwWork1.
 */
int filters_fft_spectrum(filters_t *fHdl);
/** Discrete wavelet transform. */
int filters_DWT(filters_t *fHdl);
/** Median filter with moving window size n. */
int filters_median(filters_t *fHdl, size_t n);
/** Trapezoidal filter as in Knoll NIMA 345(1994) 337-345.
 * @param[in] k rise time.
 * @param[in] l delay of peak.
 * l-k is the flat-top duration.
 * @param[in] M decay time constant (in number of samples) of the input
 * pulse.  Set M=-1.0 to deal with a step-like input function.
 */
int filters_trapezoidal(filters_t *fHdl, size_t k, size_t l, double M);

#endif // __FILTERS_H__
