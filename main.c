#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#include "common.h"
#include "filters.h"
#include "runScriptNGetConfig.h"
#include "hdf5rawWaveformIo.h"

RAW_WAVEFORM_BASE_TYPE *waveformBuf;

int main(int argc, char **argv)
{
    size_t i, j, iCh, iEvent=0, nEvents=0, iFrame=0, frameSize, nEventsInFile;
    char *inFileName;
    
    struct hdf5rawWaveformIo_waveform_file *waveformFile;
    struct waveform_attribute waveformAttr;
    struct hdf5rawWaveformIo_waveform_event waveformEvent;

    config_parameters_t *cParms;
    filters_t *fHdl;

    if(argc<2) {
        fprintf(stderr, "%s inFileName [iEvent] [nEvents]\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    inFileName = argv[1];
    waveformFile = hdf5rawWaveformIo_open_file_for_read(inFileName);
    if(argc>2)
        iEvent = atol(argv[2]);
    if(argc>3)
        nEvents = atol(argv[3]);

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

    cParms = get_config_parameters();

    fHdl = filters_init_for_convolution(frameSize, cParms->filter_respLen);
    filters_SavitzkyGolay(fHdl, cParms->filter_SavitzkyGolay_poly_order,
                          cParms->filter_SavitzkyGolay_derivative_degree);
    filters_raisedCosine(fHdl);

    iCh = 0;
    for(waveformEvent.eventId = iEvent; waveformEvent.eventId < iEvent + nEvents;
        waveformEvent.eventId++) {
        hdf5rawWaveformIo_read_event(waveformFile, &waveformEvent);
        for(iFrame = 0; iFrame < waveformAttr.nFrames; iFrame++) {
            for(i=0; i<frameSize; i++) {
                fHdl->inWav[i] = 101.5 - waveformBuf[iCh * waveformFile->nPt + iFrame * frameSize + i];
            }
            filters_convolute(fHdl);
            for(i=0; i<frameSize; i++) {
                printf("%24.16e %24.16e\n", fHdl->inWav[i], fHdl->outWav[i]);
            }
            printf("\n");
        }
        printf("\n");
    }

/*
    for(waveformEvent.eventId = iEvent; waveformEvent.eventId < iEvent + nEvents;
        waveformEvent.eventId++) {
        hdf5rawWaveformIo_read_event(waveformFile, &waveformEvent);

        for(i = 0; i < waveformFile->nPt; i++) {
            printf("%24.16e ", waveformAttr.dt*(i%frameSize));
            j = 0;
            for(iCh=0; iCh<SCOPE_NCH; iCh++) {
                if((1<<iCh) & waveformAttr.chMask) {
                    printf("%24.16e ", (waveformBuf[j * waveformFile->nPt + i]
                                        - waveformAttr.yoff[iCh]) * waveformAttr.ymult[iCh]);
                    j++;
                }
            }
            printf("\n");
            if((i+1) % frameSize == 0)
                printf("\n");
        }
        printf("\n");
    }
*/
    free(waveformBuf);
    hdf5rawWaveformIo_close_file(waveformFile);
    filters_close(fHdl);
    free_config_parameters(cParms);

    return EXIT_SUCCESS;
}
