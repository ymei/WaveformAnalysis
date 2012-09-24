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
#include "hdf5rawWaveformIo.h"

RAW_WAVEFORM_BASE_TYPE *waveformBuf;

int main(int argc, char **argv)
{
    size_t i, j, iCh, iEvent=0, nEvents=0, iFrame=0, frameSize, nEventsInFile;
    size_t nWaves2print=0, iWaves2print;
    char *configFileName, *inFileName;
    
    struct hdf5rawWaveformIo_waveform_file *waveformFile;
    struct waveform_attribute waveformAttr;
    struct hdf5rawWaveformIo_waveform_event waveformEvent;

    config_parameters_t *cParms;
    peakfinder_t *pfHdl;
    ANALYSIS_WAVEFORM_BASE_TYPE *inWav;

    if(argc<4) {
        fprintf(stderr, "%s config.scm inFileName iCh [iEvent] [nEvents] [nWaves2print]\n",argv[0]);
        return EXIT_FAILURE;
    }

    configFileName = argv[1];
    inFileName = argv[2];
    waveformFile = hdf5rawWaveformIo_open_file_for_read(inFileName);
    iCh = atol(argv[3]);
    if(argc>4)
        iEvent = atol(argv[4]);
    if(argc>4)
        nEvents = atol(argv[5]);
    if(argc>5)
        nWaves2print = atol(argv[6]);

    hdf5rawWaveformIo_read_waveform_attribute_in_file_header(waveformFile, &waveformAttr);
    fprintf(stderr, "waveform_attribute:\n"
            "     chMask  = 0x%02x\n"
            "     nPt     = %zd\n"
            "     nFrames = %zd\n"
            "     dt      = %g\n"
            "     t0      = %g\n"
            "     ymult   = %g %g %g %g\n"
            "     yoff    = %g %g %g %g\n"
            "     yzero   = %g %g %g %g\n",
            waveformAttr.chMask, waveformAttr.nPt, waveformAttr.nFrames, waveformAttr.dt,
            waveformAttr.t0, waveformAttr.ymult[0], waveformAttr.ymult[1], waveformAttr.ymult[2],
            waveformAttr.ymult[3], waveformAttr.yoff[0], waveformAttr.yoff[1],
            waveformAttr.yoff[2], waveformAttr.yoff[3], waveformAttr.yzero[0],
            waveformAttr.yzero[1], waveformAttr.yzero[2], waveformAttr.yzero[3]);

    nEventsInFile = hdf5rawWaveformIo_get_number_of_events(waveformFile);
    fprintf(stderr, "Number of events in file: %zd\n", nEventsInFile);
    if(nEvents <= 0 || nEvents > nEventsInFile) nEvents = nEventsInFile;
    if(waveformAttr.nFrames > 0) {
        frameSize = waveformAttr.nPt / waveformAttr.nFrames;
        fprintf(stderr, "Frame size: %zd\n", frameSize);
    } else {
        frameSize = waveformAttr.nPt;
    }

    waveformBuf = (RAW_WAVEFORM_BASE_TYPE*)malloc(waveformFile->nPt * waveformFile->nCh 
                                                  * sizeof(RAW_WAVEFORM_BASE_TYPE));
    waveformEvent.wavBuf = waveformBuf;

    inWav = (ANALYSIS_WAVEFORM_BASE_TYPE *)calloc(frameSize, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));
    cParms = get_config_parameters(configFileName);
    pfHdl = peakfinder_init(frameSize, 10, cParms);

    iWaves2print = 0;
    for(waveformEvent.eventId = iEvent; waveformEvent.eventId < iEvent + nEvents;
        waveformEvent.eventId++) {
        hdf5rawWaveformIo_read_event(waveformFile, &waveformEvent);
        for(iFrame = 0; iFrame < MAX(waveformAttr.nFrames, 1); iFrame++) {
            for(i=0; i<frameSize; i++) {
                inWav[i] = (waveformBuf[iCh * waveformFile->nPt + iFrame * frameSize + i]
                            - waveformAttr.yoff[iCh]) * waveformAttr.ymult[iCh];
            }

            peakfinder_baseline(pfHdl, inWav, 1);
            peakfinder_find(pfHdl);

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
                for(i=0; i<frameSize; i++) {
                    printf("%24.16e %24.16e %24.16e\n", waveformAttr.dt*i, pfHdl->blsWav[i],
                           pfHdl->fHdl->outWav[i]);
                }
                printf("\n\n");
                iWaves2print++;
            }

            if(nWaves2print == 0) {
                for(j=0; j<pfHdl->nPeaks; j++) {
                    printf("%zd %zd %g %g %g %g\n", waveformEvent.eventId, iFrame,
                           pfHdl->pParms[j].pTime, pfHdl->pParms[j].pIntegral,
                           pfHdl->pParms[j].pHeight, pfHdl->pParms[j].pWidth);
                }
            }
        }
    }

    free(waveformBuf);
    hdf5rawWaveformIo_close_file(waveformFile);
    free(inWav);
    peakfinder_close(pfHdl);
    free_config_parameters(cParms);

    return EXIT_SUCCESS;
}
