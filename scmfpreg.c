/*
 * (Tiny)SCHEME foreign pointer registry
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>

#include "scheme-private.h"
#include "scheme.h"
#include "common.h"
#include "scmfpreg.h"

struct scheme_fp 
{
    size_t n;
    void **list;
};

void *scheme_fp_init(scheme *sc)
{
    struct scheme_fp *reg;
    reg = calloc(1, sizeof(struct scheme_fp));
    if(reg) {
        if(sc->ext_data) { debug_printf("sc->ext_data seems to have valid contents.\n"); }
        scheme_set_external_data(sc, reg);
    } else {
        error_printf("malloc scheme_fp_reg failed.\n");
        return NULL;
    }
    reg->n = SCHEME_FP_NENTRIES;
    reg->list = (void**)calloc(1, reg->n * sizeof(void*));
    if(reg->list == NULL) {
        error_printf("malloc scheme_fp_reg_list failed.\n");
        free(reg);
        return NULL;
    }
    return reg;
}

int scheme_fp_close(scheme *sc)
{
    struct scheme_fp *reg;
    ssize_t i;
    int n=0;
    
    reg = (struct scheme_fp *)sc->ext_data;
    for(i=0; i<reg->n; i++) {
        if(*(reg->list + i))
            n++;
    }
    if(n) {
        error_printf("%d elements in the fp_list still hold data.\n", n);
        error_printf("Freeing anyway.\n");
    }
    if(reg->list) free(reg->list);
    if(reg) free(reg);
    return n;
}

scmfpid_t scheme_fp_register(scheme *sc, const void *p)
{
    struct scheme_fp *reg;
    ssize_t i, j;
    scmfpid_t ret=0;
    void **newlist;
    
    reg = (struct scheme_fp *)sc->ext_data;
    for(i=0; i<reg->n; i++) {
        if(*(reg->list + i) == NULL) {
            *(reg->list + i) = (void*)p;
            ret = (scmfpid_t)i;
            break;
        }
    }
    if(i == reg->n) { /* all slots in the list are occupied */
        newlist = (void**)realloc(reg->list, reg->n * 2 * sizeof(void*));
        if(newlist == NULL) {
            error_printf("realloc scheme_fp_reg_list failed.\n");
            ret = (scmfpid_t)(i+1);
        } else {
            j = reg->n;
            reg->n *= 2;
            reg->list = newlist;
            
            *(reg->list + j) = (void*)p;
            ret = (scmfpid_t)j;

            for(i=j+1; i<reg->n; i++) {
                *(reg->list + i) = NULL;
            }
        }   
    }   
    return ret; /* ret > reg->n signals failure */
}

scmfpid_t scheme_fp_remove(scheme *sc, scmfpid_t id)
{
    struct scheme_fp *reg;
    ssize_t i;
    scmfpid_t ret=0;
    
    reg = (struct scheme_fp *)sc->ext_data;
    i = id;
    if(i<reg->n && i>=0 && *(reg->list+i)) {
        *(reg->list+i) = NULL;
        ret = id;
    } else {
        error_printf("id = %zd doesn't contain registered pointer.\n", i);
        ret = (scmfpid_t)reg->n;
    }
    return ret; /* ret == reg->n signals failure */
}

void *scheme_fp_getp(scheme *sc, scmfpid_t id)
{
    struct scheme_fp *reg;
    ssize_t i;
    
    reg = (struct scheme_fp *)sc->ext_data;
    i = id;
    if(i<reg->n && i>=0 && *(reg->list+i)) {
        return *(reg->list+i);
    } else {
        error_printf("id = %zd doesn't contain registered pointer.\n", i);
        return NULL;
    }
    return NULL;
}

scmfpid_t scheme_fp_findp(scheme *sc, const void *p)
{
    struct scheme_fp *reg;
    ssize_t i, j;
    
    reg = (struct scheme_fp *)sc->ext_data;
    for(i=0; i<reg->n; i++) {
        if(*(reg->list + i) == p) {
            j = i;
        }
    }
    if(i == reg->n) {
        j = i;
        error_printf("pointer %p wasn't found in the registry.\n", p);
    }
    return (scmfpid_t)j; /* ret == reg->n signals error */
}

#ifdef SCMFPREG_DEBUG_ENABLEMAIN
int main(int argc, char **argv)
{
    scheme *sc;
    ssize_t i;

    /* init the interpreter */
    if ((sc = scheme_init_new()) == NULL) {
        error_printf("Could not initialize TinyScheme!");
        return EXIT_FAILURE;
    }
    sc->ext_data = NULL;

    scheme_set_input_port_file(sc, stdin);
    scheme_set_output_port_file(sc, stderr);

    printf("reg = %p\n", scheme_fp_init(sc));
    printf("reg id = %d\n", scheme_fp_register(sc, &i));
    printf("%p at reg id = %d\n", scheme_fp_getp(sc, 0), 0);
    printf("remove reg id = %d, ret = %d\n", 0, scheme_fp_remove(sc, 0));

    scheme_load_named_file(sc,stdin,0);
    scheme_fp_close(sc);
    
    scheme_deinit(sc);
    return EXIT_SUCCESS;
}
#endif
