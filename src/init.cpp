/* most of the content here is genereted by calling
 * `Rcpp::compileAttributes()` */

#define INCLUDE_RCPP
#include "tmb_includes.h"

/* FROM: COMPILEATTRIBUTES */
using namespace Rcpp;

// get_VA_funcs
SEXP get_VA_funcs(Rcpp::List data, Rcpp::List parameters);
RcppExport SEXP _survTMB_get_VA_funcs(SEXP dataSEXP, SEXP parametersSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
  Rcpp::traits::input_parameter< Rcpp::List >::type parameters(parametersSEXP);
  rcpp_result_gen = Rcpp::wrap(get_VA_funcs(data, parameters));
  return rcpp_result_gen;
  END_RCPP
}
// VA_funcs_eval_lb
double VA_funcs_eval_lb(SEXP p, SEXP par);
RcppExport SEXP _survTMB_VA_funcs_eval_lb(SEXP pSEXP, SEXP parSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  rcpp_result_gen = Rcpp::wrap(VA_funcs_eval_lb(p, par));
  return rcpp_result_gen;
  END_RCPP
}
// VA_funcs_eval_grad
Rcpp::NumericVector VA_funcs_eval_grad(SEXP p, SEXP par);
RcppExport SEXP _survTMB_VA_funcs_eval_grad(SEXP pSEXP, SEXP parSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  rcpp_result_gen = Rcpp::wrap(VA_funcs_eval_grad(p, par));
  return rcpp_result_gen;
  END_RCPP
}
// VA_funcs_eval_hess
Rcpp::NumericMatrix VA_funcs_eval_hess(SEXP p, SEXP par);
RcppExport SEXP _survTMB_VA_funcs_eval_hess(SEXP pSEXP, SEXP parSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  rcpp_result_gen = Rcpp::wrap(VA_funcs_eval_hess(p, par));
  return rcpp_result_gen;
  END_RCPP
}
// VA_funcs_eval_hess_sparse
Rcpp::List VA_funcs_eval_hess_sparse(SEXP p, SEXP par);
RcppExport SEXP _survTMB_VA_funcs_eval_hess_sparse(SEXP pSEXP, SEXP parSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  rcpp_result_gen = Rcpp::wrap(VA_funcs_eval_hess_sparse(p, par));
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP run_testthat_tests();

static const R_CallMethodDef CallEntries[] = {
  /* END: COMPILEATTRIBUTES */
  TMB_CALLDEFS,
  /* FROM: COMPILEATTRIBUTES */
  {"_survTMB_get_VA_funcs", (DL_FUNC) &_survTMB_get_VA_funcs, 2},
  {"_survTMB_VA_funcs_eval_lb", (DL_FUNC) &_survTMB_VA_funcs_eval_lb, 2},
  {"_survTMB_VA_funcs_eval_grad", (DL_FUNC) &_survTMB_VA_funcs_eval_grad, 2},
  {"_survTMB_VA_funcs_eval_hess", (DL_FUNC) &_survTMB_VA_funcs_eval_hess, 2},
  {"_survTMB_VA_funcs_eval_hess_sparse", (DL_FUNC) &_survTMB_VA_funcs_eval_hess_sparse, 2},
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
