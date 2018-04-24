#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(fill1d)(double *vec, int *DIMP, double *parms, int *II);
extern void F77_NAME(initmod)(void (*steadyparms)(int *, double *));
extern void F77_NAME(runmod)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(runmodbw)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);

static const R_FortranMethodDef FortranEntries[] = {
  {"fill1d", (DL_FUNC) &F77_NAME(fill1d),  4},
  {"initmod", (DL_FUNC) &F77_NAME(initmod),  1},
  {"runmod", (DL_FUNC) &F77_NAME(runmod),  6},
  {"runmodbw", (DL_FUNC) &F77_NAME(runmodbw),  6},
  {NULL, NULL, 0}
};

void R_init_DDD(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
