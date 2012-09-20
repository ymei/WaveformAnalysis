#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_wavelet.h>
#include <fftw3.h>
#include "common.h"
#include "filters.h"

filters_t *filters_init(size_t n) /* input the waveform length.  The output is also of length n */
{
    filters_t *fHdl;
    fHdl = (filters_t*)malloc(sizeof(filters_t));

    fHdl->wavLen = n;
    fHdl->respLen = 0;
    fHdl->fftUsed = 0;
    
    fHdl->inWav = (ANALYSIS_WAVEFORM_BASE_TYPE*)
        calloc(fHdl->wavLen, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    fHdl->outWav = (ANALYSIS_WAVEFORM_BASE_TYPE*)
        calloc(fHdl->wavLen, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    fHdl->waveletWav = (WAVELET_BASE_TYPE*)
        calloc(fHdl->wavLen, sizeof(WAVELET_BASE_TYPE));

    fHdl->gslDWT = gsl_wavelet_alloc(gsl_wavelet_daubechies_centered, 10);
    fHdl->gslDWTWork = gsl_wavelet_workspace_alloc(fHdl->wavLen);

    return fHdl;
}

filters_t *filters_init_for_convolution(size_t n, size_t np) /* for convolution */
{
    filters_t *fHdl;

    if(np % 2 == 0) {
        error_printf("%s(): np = %zd is not odd!\n", __FUNCTION__, np);
        return NULL;
    }

    fHdl = filters_init(n);
    fHdl->respLen = np;
    fHdl->fftUsed = 1;
    fHdl->respWav = (ANALYSIS_WAVEFORM_BASE_TYPE*)
        calloc(fHdl->respLen, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    
    //fftlen = (wavlen + resplen/2);
    fHdl->fftLen = fHdl->wavLen + (((fHdl->respLen/2)<128)?128:(fHdl->wavLen/2));

    fHdl->fftWork = (FFT_BASE_TYPE*) FFTW(malloc)(sizeof(FFT_BASE_TYPE) * fHdl->fftLen);
    fHdl->fftWork1 = (FFT_BASE_TYPE*) FFTW(malloc)(sizeof(FFT_BASE_TYPE) * fHdl->fftLen);

    fHdl->fftwPlan = FFTW(plan_r2r_1d)(fHdl->fftLen, fHdl->fftWork, fHdl->fftWork,
                                       FFTW_R2HC, FFTW_MEASURE);
    fHdl->fftwPlan1 = FFTW(plan_r2r_1d)(fHdl->fftLen, fHdl->fftWork1, fHdl->fftWork1,
                                        FFTW_R2HC, FFTW_MEASURE);
    fHdl->fftwPlan2 = FFTW(plan_r2r_1d)(fHdl->fftLen, fHdl->fftWork, fHdl->fftWork,
                                        FFTW_HC2R, FFTW_MEASURE);

    return fHdl;
}

int filters_close(filters_t *fHdl)
{
    if(fHdl->outWav)
        free(fHdl->outWav);
    if(fHdl->inWav)
        free(fHdl->inWav);
    if(fHdl->fftUsed) {
        if(fHdl->respWav)
            free(fHdl->respWav);
        FFTW(destroy_plan)(fHdl->fftwPlan);
        FFTW(destroy_plan)(fHdl->fftwPlan1);
        FFTW(destroy_plan)(fHdl->fftwPlan2);
        FFTW(free)(fHdl->fftWork);
        FFTW(free)(fHdl->fftWork1);
     }
    gsl_wavelet_free(fHdl->gslDWT);
    gsl_wavelet_workspace_free(fHdl->gslDWTWork);
    if(fHdl->waveletWav)
        free(fHdl->waveletWav);
    return 0;
}

/* m: order of polynomial, np: number of points, ld: degree of derivative*/
int filters_SavitzkyGolay(filters_t *fHdl, int m, int ld)
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
{
    size_t i;
    ANALYSIS_WAVEFORM_BASE_TYPE x;
    
    for(i=0; i<fHdl->respLen; i++) {
        x = 2.0*M_PI/(ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->respLen * (ANALYSIS_WAVEFORM_BASE_TYPE)i;
        fHdl->respWav[i] = (1.0 + cos(x))/(ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->respLen;
    }
    return 0;
}

int filters_convolute(filters_t *fHdl)
{
    size_t i;
    ANALYSIS_WAVEFORM_BASE_TYPE re, im;
    
    for(i=0; i<fHdl->wavLen; i++) {
        fHdl->fftWork[i] = fHdl->inWav[i];
    }
    for(i=fHdl->wavLen; i<fHdl->fftLen; i++) {
        fHdl->fftWork[i] = 0.0;
    }

    // fill in with respwav in wrap-around order
    fHdl->fftWork1[0] = fHdl->respWav[0];
    for(i=1; i<(fHdl->respLen+1)/2; i++) {
        fHdl->fftWork1[i] = fHdl->respWav[i];
        fHdl->fftWork1[fHdl->fftLen-i] = fHdl->respWav[fHdl->respLen-i];
    }
    for(i=(fHdl->respLen+1)/2; i<(fHdl->fftLen-(fHdl->respLen+1)/2); i++) {
        fHdl->fftWork1[i] = 0.0;
    }
    
    // do fft
    FFTW(execute)(fHdl->fftwPlan);
    FFTW(execute)(fHdl->fftwPlan1);

    // multiply in complex fourier space, half-complex format
    fHdl->fftWork[0] = fHdl->fftWork[0] * fHdl->fftWork1[0] 
        / (ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->fftLen;
    for(i=1; i<fHdl->fftLen/2; i++) {
        re = fHdl->fftWork[i] * fHdl->fftWork1[i]
            - fHdl->fftWork[fHdl->fftLen-i] * fHdl->fftWork1[fHdl->fftLen-i];
        im = fHdl->fftWork[i] * fHdl->fftWork1[fHdl->fftLen-i]
            + fHdl->fftWork[fHdl->fftLen-i] * fHdl->fftWork1[i];
        fHdl->fftWork[i] = re / (ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->fftLen;
        fHdl->fftWork[fHdl->fftLen-i] = im / (ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->fftLen;
    }
    fHdl->fftWork[fHdl->fftLen/2] = fHdl->fftWork[fHdl->fftLen/2] * fHdl->fftWork1[fHdl->fftLen/2]
        / (ANALYSIS_WAVEFORM_BASE_TYPE)fHdl->fftLen;

    // ifft
    FFTW(execute)(fHdl->fftwPlan2);

    // copy the output to outwav
    memcpy(fHdl->outWav, fHdl->fftWork, fHdl->wavLen * sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    return 0;
}

int filters_DWT(filters_t *fHdl) /* discrete wavelet transform */
{
    gsl_wavelet_transform_forward(fHdl->gslDWT, fHdl->waveletWav, 1, fHdl->wavLen,
                                  fHdl->gslDWTWork);
    return 0;
}

#ifdef FILTERS_DEBUG_ENABLEMAIN
int main(int argc, char **argv)
{
    #define PLEN 256
    ANALYSIS_WAVEFORM_BASE_TYPE pulse[PLEN] = {0.0};
    size_t i;
    filters_t *fHdl;

    i = 2;
    pulse[i++]=0.0; pulse[i++]=1.0; pulse[i++]=10.0; pulse[i++]=8.0; pulse[i++]=6.0;
    pulse[i++]=4.0; pulse[i++]=2.0; pulse[i++]=1.0; pulse[i++]=0.5; pulse[i++]=0.2;
    i = 100;
    pulse[i--]=0.0; pulse[i--]=1.0; pulse[i--]=10.0; pulse[i--]=8.0; pulse[i--]=6.0;
    pulse[i--]=4.0; pulse[i--]=2.0; pulse[i--]=1.0; pulse[i--]=0.5; pulse[i--]=0.2;

    fHdl = filters_init_for_convolution(PLEN, 31);
    memcpy(fHdl->inWav, pulse, PLEN * sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    
    filters_raisedCosine(fHdl);
    filters_convolute(fHdl);

    for(i=0; i<PLEN; i++) {
        printf("%g %g\n", pulse[i], fHdl->outWav[i]);
    }
    printf("\n\n");

    filters_SavitzkyGolay(fHdl, 10, 1);
    filters_convolute(fHdl);

    for(i=0; i<PLEN; i++) {
        printf("%g %g\n", pulse[i], fHdl->outWav[i]);
    }

    filters_close(fHdl);
    return EXIT_SUCCESS;
}
#endif
