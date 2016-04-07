#ifndef __SCMFPREG_H__
#define __SCMFPREG_H__
#include "scheme.h"

/* Routines for registering and accessing foreign void* pointers into
 * the sc->ext_data, so that foreign function libraries can store
 * object handles there. */

#define SCHEME_FP_NENTRIES 1024
typedef int scmfpid_t; /* must be a signed type */

void *scheme_fp_init(scheme *sc);
int scheme_fp_close(scheme *sc);
scmfpid_t scheme_fp_register(scheme *sc, const void *p);
scmfpid_t scheme_fp_remove(scheme *sc, scmfpid_t id);
void *scheme_fp_getp(scheme *sc, scmfpid_t id);
scmfpid_t scheme_fp_findp(scheme *sc, const void *p);

#endif /* __SCMFPREG_H__ */
