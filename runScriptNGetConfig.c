#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <readline/readline.h>
#include <readline/history.h>

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
/* list of object names in scheme, for rl completion */
static char **oblist_names=NULL;

static inline char *strformat4ts(char *in)
/* Convert letters to lower case and replace _ with - */
{
    char *p;
    for (p=in ; *p ; ++p) { *p = tolower(*p); if(*p == '_') *p = '-';}
    return in;
}

scheme *tinyscheme_init(const char *fName)
{
    FILE *fp=NULL;
    scheme *sc;
    ssize_t i;

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
        /* String literals assigned in scmvar_config_table may be read-only on some systems,
           so we dup them.  They are freed in _close() */
        scmvar_config_table[i].scmSymName = strformat4ts(strdup(scmvar_config_table[i].scmSymName));
        scmvar_config_table[i].symbol = sc->vptr->mk_symbol(sc, scmvar_config_table[i].scmSymName);
        switch (scmvar_config_table[i].scmVart) {
        case scmvart_integer:
            sc->vptr->scheme_define(sc, sc->global_env, scmvar_config_table[i].symbol,
                                    sc->vptr->mk_integer(sc, 0));
            /* If want to make constant:   sc->vptr->setimmutable(symbol); */
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
    ssize_t i;

    scheme_deinit(sc);
    for(i=0; i<scmvar_config_table_n; i++) {
        if(scmvar_config_table[i].scmSymName)
            free(scmvar_config_table[i].scmSymName);
    }
    return 0;
}

config_parameters_t *get_config_parameters(const char *fName)
{
    ssize_t i;
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

static void build_completion_list(scheme *sc)
{
    pointer scp;
    ssize_t i, j;
    size_t n=128;
    if((oblist_names = calloc(n, sizeof(char*)))==NULL) {
        error_printf("calloc oblist_names error.\n");
        return;
    }
    
    j = 0;
    for(i=0; i<sc->vptr->vector_length(sc->oblist); i++) {
        scp = sc->vptr->vector_elem(sc->oblist, i);
        do {
            oblist_names[j++] = sc->vptr->symname(sc->vptr->pair_car(scp));
            if(j>=n) {
                n *= 2;
                if((oblist_names = realloc(oblist_names, n*sizeof(char*)))==NULL) {
                    error_printf("realloc oblist_names error.\n");
                    return;
                }
            }
        } while((scp = sc->vptr->pair_cdr(scp)) != sc->NIL);
    }
    oblist_names[j] = NULL;
}

static char *completion_entry_gen(const char *text, int state)
{
    static size_t list_index, len;
    char *name;
    if(!state) {
        list_index = 0;
        len = strlen(text);
    }
    while((name = oblist_names[list_index])) {
        list_index++;
        if(strncmp(name, text, len)==0)
            return(strdup(name)); /* freed by rl automatically */
    }
    return((char*)NULL);
}

static char **completion_func(const char *text, int start, int end)
{
    char **matches;
    matches = (char**)NULL;
    matches = rl_completion_matches(text, &completion_entry_gen);
    return matches;
}

int tinyscheme_interactive(scheme *sc)
{
    ssize_t i;
    size_t n;
    char *input, prompt[1024], *output;

    output = (char*)calloc(1, 1024);
    scheme_set_output_port_string(sc, output, output+1024);
    /* Allow the output buffer to be re-allocated and freed by tinyscheme */
    sc->outport->_object._port->kind |= port_srfi6;

    rl_readline_name = "wa";
    build_completion_list(sc);
    rl_attempted_completion_function = completion_func;

    /* Configure readline to auto-complete paths when the tab key is hit. */
    rl_bind_key('\t', rl_complete);
    for(i=0;;i++) {
        /* Create prompt string from user name and current working directory. */
        snprintf(prompt, sizeof(prompt), "wa %zd > ", i);
        /* Display prompt and read input (n.b. input must be freed after use)... */
        input = readline(prompt);
        /* Check for EOF. */
        if(!input) {
            printf("\n");
            break;
        }
        /* Add input to history. */
        add_history(input);
        /* Wrap with (display ...) so the evaluated results are printed
        n = strlen(input);
        dinput = (char*)calloc(1, n+11);
        strcpy(dinput, "(display ");
        strcat(dinput, input);
        dinput[n+9] = ')';
        */
        /* Reset the state of, and clear, the output string */
        sc->outport->_object._port->rep.string.curr = sc->outport->_object._port->rep.string.start;
        bzero(sc->outport->_object._port->rep.string.start,
              sc->outport->_object._port->rep.string.past_the_end
              - sc->outport->_object._port->rep.string.start);
        /* Evaluate the input */
        scheme_load_string(sc, input);
        /* Print the output */
        n = strlen(sc->outport->_object._port->rep.string.start);
        if(n) {
            printf("%s", sc->outport->_object._port->rep.string.start);
            if(*(sc->outport->_object._port->rep.string.start + n - 1) != '\n')
                printf("\n");
        }
        /* free(dinput); */
        free(input);
    }
    /* Output could have been re-allocated, so don't free it */
    /* sc->outport->_object._port->rep.string.start seems to be freed by tinyscheme */
    free(oblist_names);
    return i;
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

    tinyscheme_interactive(sc);
    /* If not using readline for interactive session */
    /* scheme_load_named_file(sc,stdin,0); */

    tinyscheme_close(sc);
    return EXIT_SUCCESS;
}
#endif
