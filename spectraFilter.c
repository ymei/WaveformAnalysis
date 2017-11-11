/** \file
 * Filter and compute power spectra of signal.
 * Also handles Sigma-Delta Modulator (SDM) with MASH21
 * (3rd order) architecture.
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "filters.h"
#include "mreadarray.h"

/** Parameters settable from commandline */
typedef struct param
{
    double fs;   /**< sampling frequency */
    double hth;  /**< threshold for high (logic 1) */
    double lth;  /**< threshold for low  (logic 0) */
    double roffs;/**< filter roll-off frequency start */
    double roffe;/**< filter roll-off frequency end */
    int padding; /**< zero padding? */
    int sdm;     /**< operate on MASH21 SDM data? */
    ssize_t ch;  /**< select which channel in the file to compute, a channel occupies columns ch and (ch+1) */
} param_t;

param_t param_default = {
    .fs      = 25.0e6,
    .hth     = 0.8,
    .lth     = 0.2,
    .roffs   = 0.8e6,
    .roffe   = 1.0e6,
    .padding = 1,
    .sdm     = 1,
    .ch      = 0
};

void print_usage(const param_t *pm)
{
    printf("Usage:\n");
    printf("      -f sampling frequency [%g]\n", pm->fs);
    printf("      -h threshold for high (logic 1) [%g]\n", pm->hth);
    printf("      -l threshold for low  (logic 0) [%g]\n", pm->lth);
    printf("      -s filter roll-off frequency start [%g]\n", pm->roffs);
    printf("      -e filter roll-off frequency end   [%g]\n", pm->roffe);
    printf("      -p zero padding? [%d].  Spectrum is computed on filtered data\n", pm->padding);
    printf("         when padding is enabled.\n");
    printf("      -d operate on MASH21 SDM data? [%d].\n", pm->sdm);
    printf("      -c select which channel in the file to compute [%zd]\n", pm->ch);
    printf("         for SDM a channel occupies columns ch and (ch+1).\n");
    printf("      inFileName\n");
}

/** Coefficients */
double d0 = 1/16.0 / (1/6.0 * 0.5 * 0.5) - 1.0;
double d1 = 1/8.0  / (1/6.0 * 0.5 * 0.5);

int main(int argc, char **argv)
{
    int optC = 0;
    param_t pm;
    char *inFileName;

    mrdary_hdl *mhdl;
    filters_t *fhdl;
    ANALYSIS_WAVEFORM_BASE_TYPE *wav;
    size_t np;
    ssize_t i, j, k=0;
    double mu, df;

    memcpy(&pm, &param_default, sizeof(pm));
    /* parse switches */
    while((optC = getopt(argc, argv, "c:d:e:f:h:l:p:s:")) != -1) {
        switch(optC) {
        case 'c':
            pm.ch = strtoll(optarg, NULL, 0);
            break;
        case 'd':
            pm.sdm = strtol(optarg, NULL, 0);
            break;
        case 'e':
            pm.roffe = strtod(optarg, NULL);
            break;
        case 'f':
            pm.fs = strtod(optarg, NULL);
            break;
        case 'h':
            pm.hth = strtod(optarg, NULL);
            break;
        case 'l':
            pm.lth = strtod(optarg, NULL);
            break;
        case 'p':
            pm.padding = strtol(optarg, NULL, 0);
            break;
        case 's':
            pm.roffs = strtod(optarg, NULL);
            break;
        default:
            print_usage(&pm);
            return EXIT_FAILURE;
            break;
        }
    }
    argc -= optind;
    argv += optind;
    if(argc<1) {
        print_usage(&pm);
        return EXIT_FAILURE;
    }
    inFileName = argv[0];

    mhdl = mrdary_init_f(inFileName, 65536);
    np = mrdary_read_all(mhdl);
    wav = (ANALYSIS_WAVEFORM_BASE_TYPE *)calloc(np, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));

    if(pm.sdm) {
        /* select which channel (columns) in the file to compute */
        k = pm.ch;
        for(i=0; i<np; i++) {
            for(j=k*2; j<k*2+2; j++) {
                if(*mrdary_value_mn(mhdl, i, j) > pm.hth) {
                    *mrdary_value_mn(mhdl, i, j) = 1.0;
                } else *mrdary_value_mn(mhdl, i, j) = -1.0;
            }
        }
        for(i=3; i<np; i++) {
            wav[i-3] =
                *mrdary_value_mn(mhdl, i-1, k*2)
                + *mrdary_value_mn(mhdl, i-1, k*2)   * d0
                - *mrdary_value_mn(mhdl, i-2, k*2)   * d0 * 2.0
                + *mrdary_value_mn(mhdl, i-3, k*2)   * d0
                + *mrdary_value_mn(mhdl, i,   k*2+1) * d1
                - *mrdary_value_mn(mhdl, i-1, k*2+1) * d1 * 2.0
                + *mrdary_value_mn(mhdl, i-2, k*2+1) * d1;
            // wav[i-3] = 1.0 * sin(0.5 * i);
        }
        np -= 3;
    } else {
        for(i=0; i<np; i++) {
            wav[i] = *mrdary_value_mn(mhdl, i, pm.ch);
            // wav[i] = 1.0 * sin(0.5 * i) + 0.2;
        }
    }
    mrdary_free(mhdl);

    j = (ssize_t)(pm.roffe / (pm.fs/(double)np)) | 1LL;
    if(pm.padding == 0) i = 0;
    fhdl = filters_init_for_convolution(NULL, np, j);
    // memcpy(fhdl->inWav, wav, fhdl->wavLen * sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    /* make mean = 0 */
    mu = 0.0;
    for(i=0; i<np; i++) mu += wav[i];
    mu /= (double)np;
    for(i=0; i<np; i++) {
        fhdl->inWav[i] = wav[i] - mu;
    }
    fhdl->dt = 1.0/pm.fs; /* Sampling time interval, for automatic normalization. */
    df = 1.0/fhdl->dt / (double)fhdl->fftLen;
    filters_freqResp_raisedCosine(fhdl, pm.roffs/df, pm.roffe/df);
    filters_convolute(fhdl, 1);
    if(pm.padding) /* compute spectrum after filtering when padding is enabled. */
        memcpy(fhdl->inWav, fhdl->outWav, fhdl->wavLen * sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    filters_fft_spectrum(fhdl);
    printf("# %s %zd points fed into FFT, spectrum length is %zd\n", inFileName, np, np/2+1);
    printf("# df = %g when fs = %g\n", df, 1.0/fhdl->dt);
    printf("# t, wav, filtered, f, Vrms/sqrt(Hz), Vrms\n");

    for(i=0; i<np; i++) {
        printf("%g %24.16e %g %g %24.16e %24.16e\n", i * fhdl->dt, wav[i], fhdl->outWav[i]+mu,
               i * df, fhdl->fftwWork[i], fhdl->fftwWork1[i]);
    }

    printf("\n\n");

    filters_close(fhdl);
    return EXIT_SUCCESS;
}
