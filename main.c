#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#include "common.h"
#include "peakFinder.h"
#include "runScriptNGetConfig.h"

#ifndef LINE_MAX
#define LINE_MAX 4096
#endif

double *waveformBuf;

static double *read_data_file(char *fname, size_t *np);

static double *read_data_file(char *fname, size_t *np)
{
    FILE *fp;
    char buf[LINE_MAX];
    double *wav;
    size_t bufLen;
    size_t i;
    int pId;
    double v;

    if((fp = fopen(fname, "r")) == NULL) {
        perror(fname);
        return NULL;
    }
    
    if(*np > 0) {
        bufLen = *np;
    } else {
        bufLen = 1000;
    }
    
    wav = (double *) malloc(sizeof(double) * bufLen);
    
    i = 0;
    while(fgets(buf, sizeof(buf), fp)) {
        sscanf(buf, "%d%lf", &pId, &v);
        wav[i] = v; i++;
        if(i >= bufLen) {
            bufLen *= 2;
            wav = (double *) realloc(wav, sizeof(double) * bufLen);
        }
    }
    fclose(fp);
    *np = i;
    return wav;
}

int main(int argc, char **argv)
{
    size_t i, j, iCh, nPt, iEvent, nEvents;
    size_t nWaves2print=0, iWaves2print;
    char *configFileName, *inFileName;
    
    struct waveform_attribute waveformAttr;

    config_parameters_t *cParms;
    peakfinder_t *pfHdl;

    if(argc<4) {
        error_printf("%s config.scm inFileName iCh [iEvent] [nEvents] [nWaves2print]\n",
                     argv[0]);
        return EXIT_FAILURE;
    }

    configFileName = argv[1];
    inFileName = argv[2];

    iCh = atol(argv[3]);
    if(argc>4)
        iEvent = atol(argv[4]);
    if(argc>4)
        nEvents = atol(argv[5]);
    if(argc>5)
        nWaves2print = atol(argv[6]);

    nPt = 0;
    waveformBuf = read_data_file(inFileName, &nPt);
    if(nPt % 2) nPt--; /* ensure nPt is an even number */

    cParms = get_config_parameters(configFileName);
    error_printf("nPt = %zd\n", nPt);
    pfHdl = peakfinder_init(nPt, 10, cParms);

    for(i=0; i<nPt; i++) {
        pfHdl->blsWav[i] = waveformBuf[i] - 100000; /* - (-118865.0); */
        pfHdl->fHdl->inWav[i] = pfHdl->blsWav[i];
    }
/*
    filters_dofft(pfHdl->fHdl);
    for(i=0; i<(pfHdl->fHdl->fftLen/2+1); i++) {
        printf("%zd %g\n", i, pfHdl->fHdl->fftWork1[i]);
    }
*/
// #if 0
    peakfinder_find_with_zero_crossing(pfHdl);

    iWaves2print = 0;

    if(iWaves2print < nWaves2print) {
        printf("# id %zd bl %g blSD %g nPeaks %zd", iWaves2print, pfHdl->baseLine,
               pfHdl->baseLineSD, pfHdl->nPeaks);
        for(j=0; j<pfHdl->nPeaks; j++) {
            printf(" (( %g %g %g )( %g %g %g))", pfHdl->pParms[j].pStart, 
                   pfHdl->pParms[j].pTime, pfHdl->pParms[j].pEnd,
                   pfHdl->pParms[j].pHeight, pfHdl->pParms[j].pWidth,
                   pfHdl->pParms[j].pIntegral);
        }
        printf("\n");
        for(i=0; i<nPt; i++) {
            printf("%24.16e %24.16e %24.16e\n", waveformAttr.dt*i, pfHdl->blsWav[i],
                   pfHdl->fHdl->outWav[i]);
        }
        printf("\n\n");
        iWaves2print++;
    }

    if(nWaves2print == 0) {
        for(j=0; j<pfHdl->nPeaks; j++) {
            printf("%g %g %g %g\n",
                   pfHdl->pParms[j].pTime, pfHdl->pParms[j].pIntegral,
                   pfHdl->pParms[j].pHeight, pfHdl->pParms[j].pWidth);
        }
    }
//#endif
    free(waveformBuf);
    peakfinder_close(pfHdl);
    free_config_parameters(cParms);

    return EXIT_SUCCESS;
}
