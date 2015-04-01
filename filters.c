#include <string.h>
#include <sys/types.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_wavelet.h>
#include <fftw3.h>
#include "common.h"
#include "utils.h"
#include "filters.h"

filters_t *filters_init(ANALYSIS_WAVEFORM_BASE_TYPE *inWav, size_t n)
/* input the waveform length.  The output is also of length n.  n can
 * be odd in principle, but an even number is preferred */
{
    filters_t *fHdl;
    fHdl = (filters_t*)malloc(sizeof(filters_t));

    fHdl->wavLen = n;
    fHdl->respLen = 0;
    fHdl->malloced = 0;
    fHdl->fftUsed = 0;
    fHdl->fftwNThreads = FFTW_NTHREADS_DEFAULT;
    fHdl->fftwFlags = FFTW_FLAGS_DEFAULT;

    if(inWav == NULL) {
        fHdl->inWav = (ANALYSIS_WAVEFORM_BASE_TYPE*)
            calloc(fHdl->wavLen, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
        fHdl->malloced = 1;
    } else {
        fHdl->inWav = inWav;
    }
    
    fHdl->outWav = (ANALYSIS_WAVEFORM_BASE_TYPE*)
        calloc(fHdl->wavLen, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    fHdl->waveletWav = (WAVELET_BASE_TYPE*)
        calloc(fHdl->wavLen, sizeof(WAVELET_BASE_TYPE));

    fHdl->gslDWT = gsl_wavelet_alloc(gsl_wavelet_daubechies_centered, 10);
    fHdl->gslDWTWork = gsl_wavelet_workspace_alloc(fHdl->wavLen);

    return fHdl;
}

filters_t *filters_init_for_convolution(ANALYSIS_WAVEFORM_BASE_TYPE *inWav, size_t n, size_t np)
/* for convolution */
{
    filters_t *fHdl;

    if((np > 0) && (np % 2 == 0)) {
        error_printf("%s(): np = %zd is not odd!\n", __FUNCTION__, np);
        return NULL;
    }

    fHdl = filters_init(inWav, n);
    fHdl->respLen = np;
    fHdl->fftUsed = 1;
    fHdl->respWav = (ANALYSIS_WAVEFORM_BASE_TYPE*)
        calloc(fHdl->respLen, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    if(np > 0) {
        fHdl->fftLen = (fHdl->wavLen + fHdl->respLen+1); /* zero padding */
        if(fHdl->fftLen % 2) fHdl->fftLen++; /* ensure fHdl->fftLen is even */
    } else {
        fHdl->fftLen = fHdl->wavLen; /* for spectrum calculation */
    }

    if(fHdl->fftwNThreads > 0) {
        if(FFTW(init_threads)() == 0) {
            error_printf("fftw_init_threads error!\n");
        }
        FFTW(plan_with_nthreads)(fHdl->fftwNThreads);
    }
    
    fHdl->fftwWork = (FFT_BASE_TYPE*) FFTW(malloc)(sizeof(FFT_BASE_TYPE) * fHdl->fftLen);
    fHdl->fftwWork1 = (FFT_BASE_TYPE*) FFTW(malloc)(sizeof(FFT_BASE_TYPE) * fHdl->fftLen);
    fHdl->fftwWin = (FFT_BASE_TYPE*) FFTW(malloc)(sizeof(FFT_BASE_TYPE) * fHdl->fftLen);
    filters_hanning_window(fHdl); /* Hanning window as the default */
    fHdl->dt = 1.0;
    
    fHdl->fftwPlan = FFTW(plan_r2r_1d)(fHdl->fftLen, fHdl->fftwWork, fHdl->fftwWork,
                                       FFTW_R2HC, fHdl->fftwFlags);
    fHdl->fftwPlan1 = FFTW(plan_r2r_1d)(fHdl->fftLen, fHdl->fftwWork1, fHdl->fftwWork1,
                                        FFTW_R2HC, fHdl->fftwFlags);
    fHdl->fftwPlan2 = FFTW(plan_r2r_1d)(fHdl->fftLen, fHdl->fftwWork, fHdl->fftwWork,
                                        FFTW_HC2R, fHdl->fftwFlags);

    return fHdl;
}

int filters_close(filters_t *fHdl)
{
    if(fHdl->malloced) {
        if(fHdl->inWav)
            free(fHdl->inWav);
    }
    if(fHdl->outWav)
        free(fHdl->outWav);
    if(fHdl->fftUsed) {
        if(fHdl->respWav)
            free(fHdl->respWav);
        FFTW(destroy_plan)(fHdl->fftwPlan);
        FFTW(destroy_plan)(fHdl->fftwPlan1);
        FFTW(destroy_plan)(fHdl->fftwPlan2);
        FFTW(free)(fHdl->fftwWork);
        FFTW(free)(fHdl->fftwWork1);
        FFTW(free)(fHdl->fftwWin);
        if(fHdl->fftwNThreads > 0) {
            FFTW(cleanup_threads)();
            FFTW(cleanup)();
        }
    }
    gsl_wavelet_free(fHdl->gslDWT);
    gsl_wavelet_workspace_free(fHdl->gslDWTWork);
    if(fHdl->waveletWav)
        free(fHdl->waveletWav);
    return 0;
}

/* Filters should directly write into fHdl->respWav the real space
 * response waveform in wrapped around order */

int filters_SavitzkyGolay(filters_t *fHdl, int m, int ld)
/* m: order of polynomial, np: number of points, ld: degree of derivative*/
{
    int np;
    ANALYSIS_WAVEFORM_BASE_TYPE *c;
    int ipj, imj, mm, j, k, nl, nr;
    double fac, sum;
    gsl_permutation * p;
    gsl_vector *b;
    gsl_matrix *a;

    np = fHdl->respLen;
    if(np<1 || np<m-1 || np%2==0 || ld>m || np!=fHdl->respLen) {
        error_printf("%s(): improper arguments, returning...\n", __FUNCTION__);
        return 1;
    }
    
    c = calloc(np, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));

    p = gsl_permutation_alloc (m+1);
    b = gsl_vector_alloc(m+1);
    a = gsl_matrix_alloc(m+1, m+1);

    nl = np/2;
    nr = nl;

    for(ipj=0;ipj<=(m << 1);ipj++) {
        sum=(ipj ? 0.0 : 1.0);
        for(k=1;k<=nr;k++) sum += pow((double)(k),(double)(ipj));
        for(k=1;k<=nl;k++) sum += pow((double)(-k),(double)(ipj));
        mm=MIN(ipj,2*m-ipj);
        for(imj=-mm;imj<=mm;imj+=2) gsl_matrix_set(a,(ipj+imj)/2,(ipj-imj)/2,sum);
    }

    gsl_linalg_LU_decomp(a, p, &k);
    for (j=0;j<m+1;j++) gsl_vector_set(b,j,0.0);
    gsl_vector_set(b,ld,1.0);
    
    gsl_linalg_LU_solve (a, p, b, b);

    for(k = -nl;k<=nr;k++) {
        sum = gsl_vector_get(b,0);
        fac = 1.0;
        for (mm=1;mm<=m;mm++) sum += gsl_vector_get(b,mm)*(fac *= k);

        j=(np-k) % np;
        c[j]=sum; // c is in wraparound order, convenient for fft convolute

        // c[nl + k] = sum;
    }
    memcpy(fHdl->respWav, c, np * sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
/*
    for(j=0; j<np; j++) {
        fprintf(stderr, "%g\n", c[j]);
    }
*/    
    gsl_vector_free(b);
    gsl_matrix_free(a);
    gsl_permutation_free(p);

/*
    for(k=nl; k<wavlen - nr; k++) {
        sum = 0.0;
        for(j=0; j<np; j++) {
            sum += c[j] * inwav[k+j-nl];
            outwav[k] = sum;
        }
    }
*/
    free(c);
    return 0;
}

int filters_raisedCosine(filters_t *fHdl)
/* outermost bin, i=(fHdl->respLen-1)/2, is 0.0 */
{
    ssize_t i;
    ANALYSIS_WAVEFORM_BASE_TYPE x;
    
    for(i=0; i<(fHdl->respLen+1)/2; i++) { /* positive side */
        x = 2.0*M_PI/(ANALYSIS_WAVEFORM_BASE_TYPE)(fHdl->respLen-1)*(ANALYSIS_WAVEFORM_BASE_TYPE)i;
        fHdl->respWav[i] = (1.0 + cos(x))/(ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->respLen;
    }
    for(i=-((fHdl->respLen-1)/2); i<0; i++) { /* negative side */
        /* i=-(fHdl->respLen-1)/2 cast to wrong value */
        x = 2.0*M_PI/(ANALYSIS_WAVEFORM_BASE_TYPE)(fHdl->respLen-1)*(ANALYSIS_WAVEFORM_BASE_TYPE)i;
        fHdl->respWav[fHdl->respLen+i] = (1.0 + cos(x))/(ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->respLen;
    }
    return 0;
}

int filters_convolute(filters_t *fHdl)
{
    size_t i;
    ANALYSIS_WAVEFORM_BASE_TYPE re, im;
    
    for(i=0; i<fHdl->wavLen; i++) {
        fHdl->fftwWork[i] = fHdl->inWav[i];
    }
    for(i=fHdl->wavLen; i<fHdl->fftLen; i++) {
        fHdl->fftwWork[i] = 0.0;
    }

    // fill in with respwav in wrap-around order
    fHdl->fftwWork1[0] = fHdl->respWav[0];
    for(i=1; i<(fHdl->respLen+1)/2; i++) {
        fHdl->fftwWork1[i] = fHdl->respWav[i];
        fHdl->fftwWork1[fHdl->fftLen-i] = fHdl->respWav[fHdl->respLen-i];
    }
    for(i=(fHdl->respLen+1)/2; i<(fHdl->fftLen-(fHdl->respLen+1)/2); i++) {
        fHdl->fftwWork1[i] = 0.0;
    }
    
    // do fft
    FFTW(execute)(fHdl->fftwPlan);
    FFTW(execute)(fHdl->fftwPlan1);

    // multiply in complex fourier space, half-complex format
    fHdl->fftwWork[0] = fHdl->fftwWork[0] * fHdl->fftwWork1[0] 
        / (ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->fftLen;
    for(i=1; i<fHdl->fftLen/2; i++) {
        re = fHdl->fftwWork[i] * fHdl->fftwWork1[i]
            - fHdl->fftwWork[fHdl->fftLen-i] * fHdl->fftwWork1[fHdl->fftLen-i];
        im = fHdl->fftwWork[i] * fHdl->fftwWork1[fHdl->fftLen-i]
            + fHdl->fftwWork[fHdl->fftLen-i] * fHdl->fftwWork1[i];
        fHdl->fftwWork[i] = re / (ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->fftLen;
        fHdl->fftwWork[fHdl->fftLen-i] = im / (ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->fftLen;
    }
    fHdl->fftwWork[fHdl->fftLen/2] = fHdl->fftwWork[fHdl->fftLen/2] * fHdl->fftwWork1[fHdl->fftLen/2]
        / (ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->fftLen;

    // ifft
    FFTW(execute)(fHdl->fftwPlan2);

    // copy the output to outwav
    memcpy(fHdl->outWav, fHdl->fftwWork, fHdl->wavLen * sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    return 0;
}

int filters_hanning_window(filters_t *fHdl)
{
    size_t i;
    
    fHdl->fftwS1 = 0.0; fHdl->fftwS2 = 0.0;
    for(i=0; i<fHdl->fftLen; i++) {
        fHdl->fftwWin[i] = 0.5 * (1.0 - cos(2*M_PI*i/(double)fHdl->fftLen));
        fHdl->fftwS1 += fHdl->fftwWin[i];
        fHdl->fftwS2 += fHdl->fftwWin[i] * fHdl->fftwWin[i];
    }

    return 0;
}

int filters_fft_spectrum(filters_t *fHdl)
/* compute fft and store the power spectrum into fftwWork(s), normalized.
 * This function handles both odd and even number of input points.
 * The resulting power spectrum has the length (int)((n+1)/2).
 * Spectra density (linearized) is stored in fftwWork,
 * Spectrum (linearized) is stored in fftwWork1.
 */
{
    size_t i;

    for(i=0; i<fHdl->fftLen; i++) {
        fHdl->fftwWork[i] = fHdl->inWav[i] * fHdl->fftwWin[i];
    }
    FFTW(execute)(fHdl->fftwPlan);

    /* Compute linearized power spectrum into fftwWork1, in [V] for example, normalized.
     * Total length should be (int)((n+1)/2) */
    fHdl->fftwWork1[0] = fHdl->fftwWork[0] / fHdl->fftwS1;
    for(i=1; i<(fHdl->fftLen+1)/2; i++) {
        fHdl->fftwWork1[i] = hypot(fHdl->fftwWork[i],
                                   fHdl->fftwWork[fHdl->fftLen - i]) * sqrt(2.0) / fHdl->fftwS1;
    }
    /* For spectra density, normalization should be * sqrt(2.0 * dt / fftwS2) */
    for(i=0; i<(fHdl->fftLen+1)/2; i++) {
        fHdl->fftwWork[i] = fHdl->fftwWork1[i] * fHdl->fftwS1 * sqrt(fHdl->dt / fHdl->fftwS2);
    }
    return 0;
}

int filters_DWT(filters_t *fHdl) /* discrete wavelet transform */
{
    gsl_wavelet_transform_forward(fHdl->gslDWT, fHdl->waveletWav, 1, fHdl->wavLen,
                                  fHdl->gslDWTWork);
    return 0;
}

int filters_median(filters_t *fHdl, size_t n) /* median filter with moving window size n */
{
    size_t i, mid=n/2;

    for(i=0; i<mid;i++) /* overhang at the beginning */
        fHdl->outWav[i] = quickselect(fHdl->inWav, mid+i+1, (mid+i+1)/2);
    for(i=0; i<fHdl->wavLen-1 - mid; i++)
        fHdl->outWav[mid+i] = quickselect(fHdl->inWav + i, n, n/2);
    for(i=0; i<mid; i++) /* overhang at the end */
        fHdl->outWav[fHdl->wavLen-1 - i] = quickselect(fHdl->inWav + fHdl->wavLen - (mid-i+1),
                                                       mid-i+1, (mid-i+1)/2);
    return 0;
}

int filters_trapezoidal(filters_t *fHdl, size_t k, size_t l, double M)
/* Trapezoidal filter as in Knoll NIMA 345(1994) 337-345.  k is the
 * rise time, l is the delay of peak, l-k is the flat-top duration, M
 * is the decay time constant (in number of samples) of the input
 * pulse.  Set M=-1.0 to deal with a step-like input function.
 */
{
    double s, pp;
    ssize_t i, j, jk, jl, jkl;
    double vj, vjk, vjl, vjkl, dkl;

    s = 0.0; pp = 0.0;
    
    for(i=0; i<fHdl->wavLen; i++) {
        j=i; jk = j-k; jl = j-l; jkl = j-k-l;
        vj   = j>=0   ? fHdl->inWav[j]   : fHdl->inWav[0];
        vjk  = jk>=0  ? fHdl->inWav[jk]  : fHdl->inWav[0];
        vjl  = jl>=0  ? fHdl->inWav[jl]  : fHdl->inWav[0];
        vjkl = jkl>=0 ? fHdl->inWav[jkl] : fHdl->inWav[0];

        dkl = vj - vjk - vjl + vjkl;
        pp = pp + dkl;
        if(M>=0.0) {
            s = s + pp + dkl * M;
        } else { /* infinit decay time, so the input is a step function */
            s = s + dkl;
        }
        fHdl->outWav[i] = s / (fabs(M) * (double)k);
    }
    return 0;
}

#ifdef FILTERS_DEBUG_ENABLEMAIN
int main(int argc, char **argv)
{
    #define PLEN 1503
    ANALYSIS_WAVEFORM_BASE_TYPE pulse[PLEN] = {0.0};
    size_t i;
    filters_t *fHdl;

    i = 2;
    pulse[i++]=0.0; pulse[i++]=1.0; pulse[i++]=10.0; pulse[i++]=8.0; pulse[i++]=6.0;
    pulse[i++]=4.0; pulse[i++]=2.0; pulse[i++]=1.0; pulse[i++]=0.5; pulse[i++]=0.2;
    i = 100;
    pulse[i--]=0.0; pulse[i--]=1.0; pulse[i--]=10.0; pulse[i--]=8.0; pulse[i--]=6.0;
    pulse[i--]=4.0; pulse[i--]=2.0; pulse[i--]=1.0; pulse[i--]=0.5; pulse[i--]=0.2;

    fHdl = filters_init(pulse, PLEN);
    filters_median(fHdl, 11);
    for(i=0; i<fHdl->wavLen; i++) {
        printf("%g %g\n", fHdl->inWav[i], fHdl->outWav[i]);
    }
    printf("\n\n");
    
#if 0
    for(i=0; i<PLEN; i++) {
        pulse[i] += 10.0 * cos(2.0 * M_PI/3.0 * i);
    }

    fHdl = filters_init_for_convolution(pulse, PLEN, 0);
    filters_fft_spectrum(fHdl);
    for(i=0; i<fHdl->fftLen; i++) {
        printf("%g %g\n", fHdl->fftwWork[i], fHdl->fftwWork1[i]);
    }
    printf("\n\n");
#endif
#if 0
    fHdl = filters_init_for_convolution(NULL, PLEN, 31);
    memcpy(fHdl->inWav, pulse, PLEN * sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    
    filters_raisedCosine(fHdl);
    filters_convolute(fHdl);

    for(i=0; i<fHdl->fftLen/*PLEN*/; i++) {
        // printf("%g %g\n", pulse[i], fHdl->outWav[i]);
        printf("%g %g\n", fHdl->fftwWork[i], fHdl->fftwWork1[i]);
    }
    printf("\n\n");

    filters_SavitzkyGolay(fHdl, 5, 0);
    filters_convolute(fHdl);

    for(i=0; i<fHdl->fftLen/*PLEN*/; i++) {
        // printf("%g %g\n", pulse[i], fHdl->outWav[i]);
        printf("%g %g\n", fHdl->fftwWork[i], fHdl->fftwWork1[i]);
    }
#endif

    filters_close(fHdl);
    return EXIT_SUCCESS;
}
#endif
