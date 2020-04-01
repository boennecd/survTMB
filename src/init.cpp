/* most of the content here is genereted by calling
 * `Rcpp::compileAttributes()` */

#define INCLUDE_RCPP
#include "tmb_includes.h"

/* FROM: COMPILEATTRIBUTES */
using namespace Rcpp;

void hello_world();
RcppExport SEXP _survTMB_hello_world() {
  BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  hello_world();
  return R_NilValue;
  END_RCPP
}

RcppExport SEXP run_testthat_tests();

static const R_CallMethodDef CallEntries[] = {
  /* END: COMPILEATTRIBUTES */
  TMB_CALLDEFS,
  /* FROM: COMPILEATTRIBUTES */
  {"_survTMB_hello_world", (DL_FUNC) &_survTMB_hello_world, 0},
  {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 0},
  {NULL, NULL, 0}
};

RcppExport void R_init_survTMB(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  /* END: COMPILEATTRIBUTES */
#ifdef TMB_CCALLABLES
  TMB_CCALLABLES("survTMB");
#endif
  /* FROM: COMPILEATTRIBUTES */
  R_useDynamicSymbols(dll, FALSE);
}
/* END: COMPILEATTRIBUTES */
