// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void con(void *, void *, void *, void *, void *, void *);
extern void con3(void *, void *, void *, void *, void *, void *, void *, void *);
extern void cutSearch(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void cutSearch2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void giveVC(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void giveVrowSum(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void giveWC(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"con",         (DL_FUNC) &con,          6},
    {"con3",        (DL_FUNC) &con3,         8},
    {"cutSearch",   (DL_FUNC) &cutSearch,   19},
    {"cutSearch2",  (DL_FUNC) &cutSearch2,  19},
    {"giveVC",      (DL_FUNC) &giveVC,       9},
    {"giveVrowSum", (DL_FUNC) &giveVrowSum,  9},
    {"giveWC",      (DL_FUNC) &giveWC,       9},
    {NULL, NULL, 0}
};

void R_init_rocTree(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
