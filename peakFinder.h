#ifndef __PEAKFINDER_H__
#define __PEAKFINDER_H__
#include "common.h"
#include "filters.h"

typedef struct peakfinder_handle
{
    /* public */
    ANALYSIS_WAVEFORM_BASE_TYPE *blsWav; /* baseline subtracted wave */
    size_t nPeaks;
    peak_parameters_t *pParms;
    ANALYSIS_WAVEFORM_BASE_TYPE baseLine;
    ANALYSIS_WAVEFORM_BASE_TYPE baseLineSD; /* standard deviation */
    /* private */
    size_t wavLen;
    size_t nPeaksMax;
    config_parameters_t *cParms;
    filters_t *fHdl;
} peakfinder_t;

peakfinder_t *peakfinder_init(size_t wavLen, size_t nPeaksMax, config_parameters_t *cParms);
ANALYSIS_WAVEFORM_BASE_TYPE peakfinder_baseline(peakfinder_t *pfHdl,
                                                ANALYSIS_WAVEFORM_BASE_TYPE *inWav,
                                                int inverse);
size_t peakfinder_find(peakfinder_t *pfHdl);
size_t peakfinder_find_with_zero_crossing(peakfinder_t *pfHdl);
int peakfinder_close(peakfinder_t *pfHdl);

#endif /* __PEAKFINDER_H__ */
