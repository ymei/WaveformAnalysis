#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "scheme-private.h"
#include "scheme.h"
#include "common.h"
#include "runScriptNGetConfig.h"

#include "init.scm.h"

/* Here everything should be in sync with config_parameters_t in common.h.
   Don't forget to modify get_config_parameters() accordingly. */
#define scmvar_config_table_n 7
static scmvar_config_t scmvar_config_table[scmvar_config_table_n] = {
    {"wa-filter_respLen", scmvart_integer, NULL},
    {"wa-filter_SavitzkyGolay_poly_order", scmvart_integer, NULL},
    {"wa-filter_SavitzkyGolay_derivative_degree", scmvart_integer, NULL},
    {"wa-peakFinder_hThreshold", scmvart_real, NULL},
    {"wa-peakFinder_minSep", scmvart_real, NULL},
    {"wa-peakFinder_integralFraction", scmvart_real, NULL},
    {"wa-baseLine_nSamples", scmvart_integer, NULL}
};

static inline char *strformat4ts(char *in)
/* convert letters to lower case and replace _ with - */
{
    char *p;
    for (p=in ; *p ; ++p) { *p = tolower(*p); if(*p == '_') *p = '-';}
    return in;
}

scheme *tinyscheme_init(char *fName)
{
    FILE *fp=NULL;
    scheme *sc;
    size_t i;

    /* init the interpreter */
    if ((sc = scheme_init_new()) == NULL) {
        error_printf("Could not initialize TinyScheme!");
        return NULL;
    }

    scheme_set_input_port_file(sc, stdin);
    scheme_set_output_port_file(sc, stderr);

    if(fName) {
        if((fp = fopen(fName, "r")) == NULL) {
            error_printf("Error opening config file %s.\n", fName);
            return NULL;
        }
    }

    for(i=0; i<scmvar_config_table_n; i++) {
        scmvar_config_table[i].scmSymName = strformat4ts(strdup(scmvar_config_table[i].scmSymName));
        scmvar_config_table[i].symbol = sc->vptr->mk_symbol(sc, scmvar_config_table[i].scmSymName);
        switch (scmvar_config_table[i].scmVart) {
        case scmvart_integer:
            sc->vptr->scheme_define(sc, sc->global_env, scmvar_config_table[i].symbol,
                                    sc->vptr->mk_integer(sc, 0));
            /* if want to make constant:   sc->vptr->setimmutable(symbol); */
            break;
        case scmvart_real:
            sc->vptr->scheme_define(sc, sc->global_env, scmvar_config_table[i].symbol,
                                    sc->vptr->mk_real(sc, 0.0));
            break;
        }
    }

    scheme_load_string(sc, init_scm_string);
    if(fp) scheme_load_named_file(sc, fp, fName);
    if(fp) fclose(fp);
    return sc;
}

int tinyscheme_close(scheme *sc)
{
    size_t i;

    scheme_deinit(sc);
    for(i=0; i<scmvar_config_table_n; i++) {
        if(scmvar_config_table[i].scmSymName)
            free(scmvar_config_table[i].scmSymName);
    }
    return 0;
}

config_parameters_t *get_config_parameters(char *fName)
{
    size_t i;
    scheme *sc;
    pointer rv;
    long rv_i;
    double rv_r;
    config_parameters_t *cParms;

    sc = tinyscheme_init(fName);
    cParms = (config_parameters_t *)malloc(sizeof(config_parameters_t));

    for(i=0; i<scmvar_config_table_n; i++) {
        rv = scheme_eval(sc, scmvar_config_table[i].symbol);
        switch (scmvar_config_table[i].scmVart) {
        case scmvart_integer:
            rv_i = sc->vptr->ivalue(rv);
            break;
        case scmvart_real:
            rv_r = sc->vptr->rvalue(rv);
            break;
        }
        switch (i) {
        case 0: cParms->filter_respLen = rv_i; break;
        case 1: cParms->filter_SavitzkyGolay_poly_order = rv_i; break;
        case 2: cParms->filter_SavitzkyGolay_derivative_degree = rv_i; break;
        case 3: cParms->peakFinder_hThreshold = rv_r; break;
        case 4: cParms->peakFinder_minSep = rv_r; break;
        case 5: cParms->peakFinder_integralFraction = rv_r; break;
        case 6: cParms->baseLine_nSamples = rv_i; break;
        }
    }

    tinyscheme_close(sc);
    return cParms;
}

int free_config_parameters(config_parameters_t *cParms)
{
    free(cParms);
    return 0;
}

#ifdef RUNSCRIPTNGETCONFIG_DEBUG_ENABLEMAIN
int main(int argc, char **argv)
{
    scheme *sc;
    pointer rv;

    sc = tinyscheme_init(NULL);

    rv = scheme_eval(sc, sc->vptr->mk_symbol(sc, "wa-filter-resplen"));
    if(rv != sc->NIL) {
        printf("%ld\n", sc->vptr->ivalue(rv));
    }

    scheme_load_named_file(sc,stdin,0);

    tinyscheme_close(sc);
    return EXIT_SUCCESS;
}
#endif
