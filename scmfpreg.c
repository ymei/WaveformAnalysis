#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>

#include "scheme-private.h"
#include "scheme.h"
#include "common.h"
#include "scmfpreg.h"

/**
 * Foreign pointer registry struct to be set to sc->ext_data.
 */
struct scheme_fp 
{
    size_t n;    /**< length of the list of elements */
    void **list; /**< list of elements */
};

static pointer scheme_fp_get_value(scheme *sc, pointer args);

/** Initialize the foreign pointer registry in the existing sc */
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
    /* register helper functions */
    sc->vptr->scheme_define(sc, sc->global_env, sc->vptr->mk_symbol(sc, "fp-get-value"),
                            sc->vptr->mk_foreign_func(sc, scheme_fp_get_value));
    return reg;
}

/** Destroy the foreign pointer registry.  Should be called after
 * removing all elements from the registry, and before destroying
 * sc */
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

/** Register a new pointer p into the registry.
 * @param[in] p pointer to be registered.
 * @return the index in the list of elements. */
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

/** Remove the registration of pointer at index rid
 * @param[in] rid
 */
scmfpid_t scheme_fp_remove(scheme *sc, scmfpid_t rid)
{
    struct scheme_fp *reg;
    ssize_t i;
    scmfpid_t ret=0;
    
    reg = (struct scheme_fp *)sc->ext_data;
    i = rid;
    if(i<reg->n && i>=0 && *(reg->list+i)) {
        *(reg->list+i) = NULL;
        ret = rid;
    } else {
        error_printf("id = %zd doesn't contain registered pointer.\n", i);
        ret = (scmfpid_t)reg->n;
    }
    return ret; /* ret == reg->n signals failure */
}

/** Retrieve the pointer at index rid
 * @param[in] rid
 * @return pointer
 */
void *scheme_fp_getp(scheme *sc, scmfpid_t rid)
{
    struct scheme_fp *reg;
    ssize_t i;
    
    reg = (struct scheme_fp *)sc->ext_data;
    i = rid;
    if(i<reg->n && i>=0 && *(reg->list+i)) {
        return *(reg->list+i);
    } else {
        error_printf("id = %zd doesn't contain registered pointer.\n", i);
        return NULL;
    }
    return NULL;
}

/** Look up if pointer p is in the registry
 * @param[in] p pointer
 * @return index where p resides.  returns reg->n when p is not in the registry.
 */
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

/****************************** Data accessors ********************************/
/** C foreign function to be presented as "fp-get-value" in scheme */
static pointer scheme_fp_get_value(scheme *sc, pointer args)
{
    scmfpid_t rid;
    pointer p1, p2, p3, p4;
    size_t i;
    char *a;
    
    if(args!=sc->NIL) {
        p1 = sc->vptr->pair_car(args);
        p2 = sc->vptr->pair_cdr(args);
        if(p2!=sc->NIL) {
            p3 = sc->vptr->pair_car(p2);
            p4 = sc->vptr->pair_cdr(p2);
        } else {
            goto error;
        }
        if(sc->vptr->is_integer(p1) && sc->vptr->is_integer(p3)) {
            rid = sc->vptr->ivalue(p1);
            i = sc->vptr->ivalue(p3);
        } else {
            goto error;
        }
        a = (char*)scheme_fp_getp(sc, rid);
        
        if(a) {
            if(p4!=sc->NIL) {
                if(sc->vptr->is_integer(sc->vptr->pair_car(p4))) {
                    switch (sc->vptr->ivalue(sc->vptr->pair_car(p4))) {
                    case 0:
                        return sc->vptr->mk_string(sc, (a+i));
                        break;
                    case 8:
                        return sc->vptr->mk_integer(sc, *((int8_t*)a+i));
                        break;
                    case 18:
                        return sc->vptr->mk_character(sc, *(a+i));
                        break;
                    case 16:
                        return sc->vptr->mk_integer(sc, *((int16_t*)a+i));
                        break;
                    case 64:
                        return sc->vptr->mk_integer(sc, *((int64_t*)a+i));
                        break;
                    case 132:
                        return sc->vptr->mk_real(sc, *((float*)a+i));
                        break;
                    case 164:
                        return sc->vptr->mk_real(sc, *((double*)a+i));
                        break;
                    case 32:
                    default:
                        return sc->vptr->mk_integer(sc, *((int32_t*)a+i));
                        break;
                    }
                } else {
                    goto error;
                }
            } else { /* default int32_t */
                return sc->vptr->mk_integer(sc, *((int32_t*)a+i));
            }
        }
    }
error:
    sc->vptr->putstr(sc, "Incorrect arguments.  Should be (fp-get-value rid i type)\n"
                         "type could be omitted (int32), or 0:string, 18:char,\n{8, 16, 32, 64}:int, 132:float, 164:double\n");
    return sc->NIL;
}

#ifdef SCMFPREG_DEBUG_ENABLEMAIN
int main(int argc, char **argv)
{
    scheme *sc;
    char *a="Hello.";
    ssize_t i = 0x12345678;

    /* init the interpreter */
    if ((sc = scheme_init_new()) == NULL) {
        error_printf("Could not initialize TinyScheme!");
        return EXIT_FAILURE;
    }
    sc->ext_data = NULL;

    scheme_set_input_port_file(sc, stdin);
    scheme_set_output_port_file(sc, stderr);

    printf("reg = %p\n", scheme_fp_init(sc));
    printf("reg id = %ld\n", scheme_fp_register(sc, &i));
    scheme_fp_register(sc, a);
    printf("%p at reg id = %ld\n", scheme_fp_getp(sc, 0L), 0L);
//    printf("remove reg id = %ld, ret = %ld\n", 0L, scheme_fp_remove(sc, 0L));

    scheme_load_named_file(sc,stdin,0);
    scheme_fp_close(sc);
    
    scheme_deinit(sc);
    return EXIT_SUCCESS;
}
#endif
