/** Compute power spectra of Sigma-Delta Modulator (SDM) with MASH21
 * (3rd order) architecture.
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "filters.h"
#include "mreadarray.h"

double hth = 3.0; /**< threshold for high (logic 1) */
double lth = 1.0; /**< threshold for low  (logic 0) */

/** Coefficients */
double d0 = 1/16.0 / (1/6.0 * 0.5 * 0.5) - 1.0;
double d1 = 1/8.0  / (1/6.0 * 0.5 * 0.5);

int main(int argc, char **argv)
{
    mrdary_hdl *mhdl;
    filters_t *fhdl;
    ANALYSIS_WAVEFORM_BASE_TYPE *wav;
    size_t np;
    ssize_t i, j;
    double df;

    mhdl = mrdary_init_f(argv[1], 65536);
    np = mrdary_read_all(mhdl);
    wav = (ANALYSIS_WAVEFORM_BASE_TYPE *)calloc(np, sizeof(ANALYSIS_WAVEFORM_BASE_TYPE));

    for(i=0; i<np; i++) {
        for(j=0; j<2; j++) {
            if(*mrdary_value_mn(mhdl, i, j) > hth) {
                *mrdary_value_mn(mhdl, i, j) = 1.0;
            } else *mrdary_value_mn(mhdl, i, j) = -1.0;
        }
    }
    for(i=3; i<np; i++) {
        wav[i-3] =
              *mrdary_value_mn(mhdl, i-1, 0)
            + *mrdary_value_mn(mhdl, i-1, 0) * d0
            - *mrdary_value_mn(mhdl, i-2, 0) * d0 * 2.0
            + *mrdary_value_mn(mhdl, i-3, 0) * d0
            + *mrdary_value_mn(mhdl, i,   1) * d1
            - *mrdary_value_mn(mhdl, i-1, 1) * d1 * 2.0
            + *mrdary_value_mn(mhdl, i-2, 1) * d1;
        // wav[i-3] = 1.0 * sin(5. * i);
    }
    mrdary_free(mhdl);

    np -= 3;

    fhdl = filters_init_for_convolution(wav, np, 0);
    fhdl->dt = 1.0/25.6e6; /* fs = 25.6MHz, for automatic normalization */
    df = 1.0/fhdl->dt / 2.0 / (double)((np+1)/2);
    filters_fft_spectrum(fhdl);
    printf("# %s %zd points fed into FFT, spectrum length is %zd\n", argv[1], np, (np+1)/2);
    printf("# df = %g when fs = %g\n", df, 1.0/fhdl->dt);
    printf("# t, wav, f, Vrms/sqrt(Hz), Vrms\n");

    for(i=0; i<np; i++) {
        printf("%g %24.16e %g %24.16e %24.16e\n", i * fhdl->dt, wav[i], i * df,
               fhdl->fftwWork[i], fhdl->fftwWork1[i]);
    }

    printf("\n\n");
    
    filters_close(fhdl);
    return EXIT_SUCCESS;
}
