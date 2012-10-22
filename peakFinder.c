#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "common.h"
#include "peakFinder.h"

peakfinder_t *peakfinder_init(size_t wavLen, size_t nPeaksMax, config_parameters_t *cParms)
{
    peakfinder_t *pfHdl;

    pfHdl = (peakfinder_t *)malloc(sizeof(peakfinder_t));
    pfHdl->wavLen = wavLen;
    pfHdl->blsWav = (ANALYSIS_WAVEFORM_BASE_TYPE *)
        malloc(pfHdl->wavLen * sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    pfHdl->nPeaksMax = nPeaksMax;
    pfHdl->pParms = (peak_parameters_t *)calloc(pfHdl->nPeaksMax, sizeof(peak_parameters_t));
    pfHdl->nPeaks = 0;

    pfHdl->cParms = cParms;
    pfHdl->fHdl = filters_init_for_convolution(pfHdl->wavLen, cParms->filter_respLen);
    // filters_raisedCosine(pfHdl->fHdl);
    filters_SavitzkyGolay(pfHdl->fHdl,
                          cParms->filter_SavitzkyGolay_poly_order,
                          cParms->filter_SavitzkyGolay_derivative_degree);
    return pfHdl;
}

int peakfinder_close(peakfinder_t *pfHdl)
{
    if(pfHdl->blsWav)
        free(pfHdl->blsWav);
    if(pfHdl->pParms)
        free(pfHdl->pParms);
    filters_close(pfHdl->fHdl);
    free(pfHdl);
    return 0;
}

ANALYSIS_WAVEFORM_BASE_TYPE peakfinder_baseline(peakfinder_t *pfHdl,
                                                ANALYSIS_WAVEFORM_BASE_TYPE *inWav,
                                                int inverse)
{
    size_t i, nSamples;

    nSamples = pfHdl->cParms->baseLine_nSamples;

    pfHdl->baseLine = 0.0;
    pfHdl->baseLineSD = 0.0;
    for(i=0; i<nSamples; i++) {
        pfHdl->baseLine += inWav[i];
    }
    pfHdl->baseLine /= (ANALYSIS_WAVEFORM_BASE_TYPE)nSamples;
    for(i=0; i<nSamples; i++) {
        pfHdl->baseLineSD += (inWav[i] - pfHdl->baseLine) * (inWav[i] - pfHdl->baseLine);
    }
    pfHdl->baseLineSD = sqrt(pfHdl->baseLineSD / (nSamples - 1.0));

    for(i=0; i<pfHdl->wavLen; i++) {
        if(inverse)
            pfHdl->blsWav[i] = pfHdl->baseLine - inWav[i];
        else
            pfHdl->blsWav[i] = inWav[i] - pfHdl->baseLine;
    }

    return pfHdl->baseLine;
}

size_t peakfinder_find(peakfinder_t *pfHdl)
{
    size_t i, j, l, r;
    int p=0, inPeak;
    ANALYSIS_WAVEFORM_BASE_TYPE max;

    for(i=0; i<pfHdl->wavLen; i++)
        pfHdl->fHdl->inWav[i] = pfHdl->blsWav[i];
    filters_convolute(pfHdl->fHdl);

    /*
     * Allow a minimum seperation (minsep) between two peaks.  If two
     * peaks are within the minsep window, and the second peak is
     * higher than the previous peak, the previous peak is discarded.
     *
     * The observation is: One high peak is always followed by a
     * smaller ringing peak within an almost fixed window.
     */

    inPeak = 0;
    pfHdl->nPeaks = 0;
    for(i=0; i<pfHdl->wavLen; i++) {
        if(inPeak == 0 && pfHdl->fHdl->outWav[i] >= pfHdl->cParms->peakFinder_hThreshold) {
            /* possible beginning of a peak (upwards peak) */
            max = pfHdl->cParms->peakFinder_hThreshold;
            p = 0;
            l = i;
            inPeak = 1;
        }
        if(inPeak == 1 && pfHdl->fHdl->outWav[i] < pfHdl->cParms->peakFinder_hThreshold) {
            /* getting out of a peak */
            inPeak = 0;
            r = i;
            if(pfHdl->nPeaks > 0) { /* compare with the previous peak */
                if((p - pfHdl->pParms[pfHdl->nPeaks-1].pTime) < pfHdl->cParms->peakFinder_minSep) {
                    /* if >, peak is added */
                    if(max > pfHdl->pParms[pfHdl->nPeaks-1].pHeight) {
                        /* the previous peak should be discarded */
                        pfHdl->nPeaks--;
                    } else { /* this peak should not be added */
                        goto GOTO_NO_ADD_PEAK;
                    }
                }
            }

            if(pfHdl->nPeaks < pfHdl->nPeaksMax) {
                pfHdl->pParms[pfHdl->nPeaks].pHeight = max;
                pfHdl->pParms[pfHdl->nPeaks].pTime = p;

                for(j=l; j<r; j++) {
                    if(pfHdl->fHdl->outWav[j] > pfHdl->cParms->peakFinder_integralFraction * max) {
                        pfHdl->pParms[pfHdl->nPeaks].pStart = j;
                        break;
                    }
                }
                for(j=r; j>=l; j--) {
                    if(pfHdl->fHdl->outWav[j] > pfHdl->cParms->peakFinder_integralFraction * max) {
                        pfHdl->pParms[pfHdl->nPeaks].pEnd = j;
                        break;
                    }
                }
                pfHdl->pParms[pfHdl->nPeaks].pHeight = 0.0;
                pfHdl->pParms[pfHdl->nPeaks].pIntegral = 0.0;
                for(j=pfHdl->pParms[pfHdl->nPeaks].pStart; j<=pfHdl->pParms[pfHdl->nPeaks].pEnd;
                    j++) {
                    pfHdl->pParms[pfHdl->nPeaks].pIntegral += pfHdl->blsWav[j];
                    if(pfHdl->blsWav[j] > pfHdl->pParms[pfHdl->nPeaks].pHeight) {
                        pfHdl->pParms[pfHdl->nPeaks].pHeight = pfHdl->blsWav[j];
                        pfHdl->pParms[pfHdl->nPeaks].pTime = j;
                    }
                }
                /* FWHM */
                for(j=pfHdl->pParms[pfHdl->nPeaks].pStart; j<=pfHdl->pParms[pfHdl->nPeaks].pEnd;
                    j++) {
                    if(pfHdl->blsWav[j] >= pfHdl->pParms[pfHdl->nPeaks].pHeight/2.0) {
                        l=j;
                        break;
                    }
                }
                for(j=pfHdl->pParms[pfHdl->nPeaks].pEnd; j>=pfHdl->pParms[pfHdl->nPeaks].pStart;
                    j--) {
                    if(pfHdl->blsWav[j] >= pfHdl->pParms[pfHdl->nPeaks].pHeight/2.0) {
                        r=j;
                        break;
                    }
                }
                pfHdl->pParms[pfHdl->nPeaks].pWidth = r-l+1.0;

                pfHdl->nPeaks++;
            }
        }

    GOTO_NO_ADD_PEAK:

        if(inPeak) {
            if(pfHdl->fHdl->outWav[i] > max) {
                max = pfHdl->fHdl->outWav[i];
                p = i;
            }
        }
    }

    return pfHdl->nPeaks;
}

size_t peakfinder_find_with_zero_crossing(peakfinder_t *pfHdl)
{
    ssize_t i, j, l, r, im, im1, tl, th, tm1, tm2;
    int p=0, inPeak;
    ANALYSIS_WAVEFORM_BASE_TYPE max, max1, bl, pH, rt, ft, fwhm;

    size_t leadSkip=100;
    
    for(i=0; i<pfHdl->wavLen; i++)
        pfHdl->fHdl->inWav[i] = pfHdl->blsWav[i];
    filters_convolute(pfHdl->fHdl);

    inPeak = 0;
    pfHdl->nPeaks = 0;
    for(i=leadSkip; i<pfHdl->wavLen; i++) {
        if(inPeak == 0 && pfHdl->fHdl->outWav[i] >= pfHdl->cParms->peakFinder_hThreshold) {
            /* possible beginning of a peak (upwards peak) */
            max = pfHdl->cParms->peakFinder_hThreshold;
            p = 0;
            l = i;
            inPeak = 1;
        }
        if(inPeak == 1 && pfHdl->fHdl->outWav[i] <= 0.0) {
            /* zero-crossing */
            r = l + 3*(i-l);
            inPeak = 0;
            max1 = -DBL_MAX;
            for(j=l; j<=r; j++) {
                if(pfHdl->fHdl->outWav[j] > max) {
                    max = pfHdl->fHdl->outWav[j];
                    im = j;
                }
                if(pfHdl->blsWav[j] > max1) {
                    max1 = pfHdl->blsWav[j];
                    im1 = j;
                }
            }
            /* baseline */
            bl = 0.0;
            for(j=0; j < pfHdl->cParms->baseLine_nSamples; j++) {
                bl += pfHdl->blsWav[l];
            }
            bl /= (ANALYSIS_WAVEFORM_BASE_TYPE)(pfHdl->cParms->baseLine_nSamples);
            pH = max1 - bl;
            /* rise and fall times, FWHM*/
            tl = 0; th = 0; tm1 = 0;
            for(j=im1; j>=l; j--) {
                if(th == 0 && (pfHdl->blsWav[j] - bl < 0.9 * pH)) { th = j+1; }
                if(tl == 0 && (pfHdl->blsWav[j] - bl < 0.1 * pH)) { tl = j; }
                if(tm1 == 0 && (pfHdl->blsWav[j] - bl <= 0.5 * pH)) { tm1 = j; }
            }
            rt = th - tl;
            tl = 0; th = 0; tm2 = 0;
            for(j=im1; j<=r; j++) {
                if(th == 0 && (pfHdl->blsWav[j] - bl < 0.9 * pH)) { th = j-1; }
                if(tl == 0 && (pfHdl->blsWav[j] - bl < 0.1 * pH)) { tl = j; }
                if(tm2 == 0 && (pfHdl->blsWav[j] - bl <= 0.5 * pH)) { tm2 = j; }
            }
            ft = tl - th;
            fwhm = tm2 - tm1 + 1;

            error_printf("%zd %g %g %g %g %g\n", im1, bl, pH, rt, ft, fwhm);
        }
    }

    return pfHdl->nPeaks;
}
