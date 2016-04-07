#ifndef __COMMON_H__
#define __COMMON_H__

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define HDF5IO(name) hdf5rawWaveformIo_ ## name

#define SCOPE_NCH 16
#define SCOPE_MEM_LENGTH_MAX 12500000 /* DPO5054 default, 12.5M points maximum */
#define SCOPE_DATA_TYPE int16_t
#define SCOPE_DATA_HDF5_TYPE H5T_NATIVE_INT16

#define RAW_WAVEFORM_BASE_TYPE SCOPE_DATA_TYPE
#define ANALYSIS_WAVEFORM_BASE_TYPE double
#define FFT_BASE_TYPE double /* if this is float, FFTW should be fftwf_ */
#define FFTW(name) fftw_ ## name
#define WAVELET_BASE_TYPE double

struct waveform_attribute 
{
    uint32_t chMask;
    uint64_t nPt;     /* number of points in each event */
    uint64_t nFrames; /* number of Fast Frames in each event, 0 means off */
    double dt;
    double t0;
    double ymult[SCOPE_NCH];
    double yoff[SCOPE_NCH];
    double yzero[SCOPE_NCH];
};

typedef struct config_parameters
{
    int filter_respLen; /* response function length, should be an odd number */
    int filter_SavitzkyGolay_poly_order;
    int filter_SavitzkyGolay_derivative_degree;
    ANALYSIS_WAVEFORM_BASE_TYPE peakFinder_hThreshold; /* height threshold */
    ANALYSIS_WAVEFORM_BASE_TYPE peakFinder_minSep; /* minimum separation between peaks */
    /* fraction of peak height in the filtered waveform, in order to
     * determin start and end points to do the integration */
    ANALYSIS_WAVEFORM_BASE_TYPE peakFinder_integralFraction;
    size_t baseLine_nSamples;
} config_parameters_t;

typedef struct peak_parameters
{
    ANALYSIS_WAVEFORM_BASE_TYPE pBaseline; /* `local' baseline for the pulse */
    ANALYSIS_WAVEFORM_BASE_TYPE pHeight; /* highest point */
    ANALYSIS_WAVEFORM_BASE_TYPE pStart; /* start on the filtered peak */
    ANALYSIS_WAVEFORM_BASE_TYPE pEnd; /* end on the filtered peak */
    ANALYSIS_WAVEFORM_BASE_TYPE pTime; /* position of the highest point on the unfiltered peak */
    ANALYSIS_WAVEFORM_BASE_TYPE pWidth; /* FWHM of the unfiltered peak */
    ANALYSIS_WAVEFORM_BASE_TYPE pIntegral; /* integral of the unfiltered peak */
} peak_parameters_t;

/* utilities */
#define bitsof(x) (8*sizeof(x))

#ifdef DEBUG
  #define debug_printf(...) do { fprintf(stderr, __VA_ARGS__); fflush(stderr); \
                               } while (0)
#else
  #define debug_printf(...) ((void)0)
#endif
#define error_printf(...) do { fprintf(stderr, __VA_ARGS__); fflush(stderr); \
                             } while(0)

#define MIN(x,y) (((x)>(y))?(y):(x))
#define MAX(x,y) (((x)<(y))?(y):(x))

#ifndef strlcpy
#define strlcpy(a, b, c) do {                   \
        strncpy(a, b, (c)-1);                   \
        (a)[(c)-1] = '\0';                      \
    } while (0)
#endif

#endif /* __COMMON_H__ */
