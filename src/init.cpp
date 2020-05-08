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
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  rcpp_result_gen = Rcpp::wrap(VA_funcs_eval_hess_sparse(p, par));
  return rcpp_result_gen;
  END_RCPP
}
// get_gl_rule
Rcpp::List get_gl_rule(unsigned const n);
RcppExport SEXP _survTMB_get_gl_rule(SEXP nSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< unsigned const >::type n(nSEXP);
  rcpp_result_gen = Rcpp::wrap(get_gl_rule(n));
  return rcpp_result_gen;
  END_RCPP
}
// joint_start_baseline
arma::vec joint_start_baseline(arma::vec const& Y, arma::vec const& tstart, arma::vec const& tstop, arma::vec const& omega, arma::vec const& offsets, unsigned const n_nodes, arma::vec const& bound_knots, arma::vec const& inter_knots, bool const grad);
RcppExport SEXP _survTMB_joint_start_baseline(SEXP YSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP omegaSEXP, SEXP offsetsSEXP, SEXP n_nodesSEXP, SEXP bound_knotsSEXP, SEXP inter_knotsSEXP, SEXP gradSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< arma::vec const& >::type Y(YSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type tstart(tstartSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type tstop(tstopSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type omega(omegaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type offsets(offsetsSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type n_nodes(n_nodesSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type bound_knots(bound_knotsSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type inter_knots(inter_knotsSEXP);
  Rcpp::traits::input_parameter< bool const >::type grad(gradSEXP);
  rcpp_result_gen = Rcpp::wrap(joint_start_baseline(Y, tstart, tstop, omega, offsets, n_nodes, bound_knots, inter_knots, grad));
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
  {"_survTMB_get_gl_rule", (DL_FUNC) &_survTMB_get_gl_rule, 1},
  {"_survTMB_joint_start_baseline", (DL_FUNC) &_survTMB_joint_start_baseline, 9},
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
