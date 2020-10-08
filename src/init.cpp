/* most of the content here is genereted by calling
 * `Rcpp::compileAttributes()` */

#define INCLUDE_RCPP
#include "tmb_includes.h"

/* FROM: COMPILEATTRIBUTES */
using namespace Rcpp;

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// psqn_get_mgsm_funcs
SEXP psqn_get_mgsm_funcs(Rcpp::List data, double const eps, double const kappa, arma::vec const& b, arma::vec const& theta, arma::vec const& theta_va, int const n_nodes, std::string const& link, unsigned const max_threads);
RcppExport SEXP _survTMB_psqn_get_mgsm_funcs(SEXP dataSEXP, SEXP epsSEXP, SEXP kappaSEXP, SEXP bSEXP, SEXP thetaSEXP, SEXP theta_vaSEXP, SEXP n_nodesSEXP, SEXP linkSEXP, SEXP max_threadsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
  Rcpp::traits::input_parameter< double const >::type eps(epsSEXP);
  Rcpp::traits::input_parameter< double const >::type kappa(kappaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type theta(thetaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type theta_va(theta_vaSEXP);
  Rcpp::traits::input_parameter< int const >::type n_nodes(n_nodesSEXP);
  Rcpp::traits::input_parameter< std::string const& >::type link(linkSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type max_threads(max_threadsSEXP);
  rcpp_result_gen = Rcpp::wrap(psqn_get_mgsm_funcs(data, eps, kappa, b, theta, theta_va, n_nodes, link, max_threads));
  return rcpp_result_gen;
  END_RCPP
}
// psqn_optim_mgsm_private
Rcpp::NumericVector psqn_optim_mgsm_private(Rcpp::NumericVector val, SEXP ptr, double const rel_eps, unsigned const max_it, unsigned const n_threads, double const c1, double const c2);
RcppExport SEXP _survTMB_psqn_optim_mgsm_private(SEXP valSEXP, SEXP ptrSEXP, SEXP rel_epsSEXP, SEXP max_itSEXP, SEXP n_threadsSEXP, SEXP c1SEXP, SEXP c2SEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
  Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
  Rcpp::traits::input_parameter< double const >::type rel_eps(rel_epsSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type max_it(max_itSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type n_threads(n_threadsSEXP);
  Rcpp::traits::input_parameter< double const >::type c1(c1SEXP);
  Rcpp::traits::input_parameter< double const >::type c2(c2SEXP);
  rcpp_result_gen = Rcpp::wrap(psqn_optim_mgsm_private(val, ptr, rel_eps, max_it, n_threads, c1, c2));
  return rcpp_result_gen;
  END_RCPP
}
// psqn_optim_mgsm
Rcpp::List psqn_optim_mgsm(Rcpp::NumericVector val, SEXP ptr, double const rel_eps, unsigned const max_it, unsigned const n_threads, double const c1, double const c2, bool const use_bfgs, int const trace, double const cg_tol, bool const strong_wolfe);
RcppExport SEXP _survTMB_psqn_optim_mgsm(SEXP valSEXP, SEXP ptrSEXP, SEXP rel_epsSEXP, SEXP max_itSEXP, SEXP n_threadsSEXP, SEXP c1SEXP, SEXP c2SEXP, SEXP use_bfgsSEXP, SEXP traceSEXP, SEXP cg_tolSEXP, SEXP strong_wolfeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type val(valSEXP);
  Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
  Rcpp::traits::input_parameter< double const >::type rel_eps(rel_epsSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type max_it(max_itSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type n_threads(n_threadsSEXP);
  Rcpp::traits::input_parameter< double const >::type c1(c1SEXP);
  Rcpp::traits::input_parameter< double const >::type c2(c2SEXP);
  Rcpp::traits::input_parameter< bool const >::type use_bfgs(use_bfgsSEXP);
  Rcpp::traits::input_parameter< int const >::type trace(traceSEXP);
  Rcpp::traits::input_parameter< double const >::type cg_tol(cg_tolSEXP);
  Rcpp::traits::input_parameter< bool const >::type strong_wolfe(strong_wolfeSEXP);
  rcpp_result_gen = Rcpp::wrap(psqn_optim_mgsm(val, ptr, rel_eps, max_it, n_threads, c1, c2, use_bfgs, trace, cg_tol, strong_wolfe));
  return rcpp_result_gen;
  END_RCPP
}
// eval_psqn_mgsm
double eval_psqn_mgsm(NumericVector val, SEXP ptr, unsigned const n_threads);
RcppExport SEXP _survTMB_eval_psqn_mgsm(SEXP valSEXP, SEXP ptrSEXP, SEXP n_threadsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< NumericVector >::type val(valSEXP);
  Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type n_threads(n_threadsSEXP);
  rcpp_result_gen = Rcpp::wrap(eval_psqn_mgsm(val, ptr, n_threads));
  return rcpp_result_gen;
  END_RCPP
}
// grad_psqn_mgsm
NumericVector grad_psqn_mgsm(NumericVector val, SEXP ptr, unsigned const n_threads);
RcppExport SEXP _survTMB_grad_psqn_mgsm(SEXP valSEXP, SEXP ptrSEXP, SEXP n_threadsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< NumericVector >::type val(valSEXP);
  Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type n_threads(n_threadsSEXP);
  rcpp_result_gen = Rcpp::wrap(grad_psqn_mgsm(val, ptr, n_threads));
  return rcpp_result_gen;
  END_RCPP
}
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
// get_commutation
Rcpp::NumericMatrix get_commutation(unsigned const n, unsigned const m);
RcppExport SEXP _survTMB_get_commutation(SEXP nSEXP, SEXP mSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< unsigned const >::type n(nSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type m(mSEXP);
  rcpp_result_gen = Rcpp::wrap(get_commutation(n, m));
  return rcpp_result_gen;
  END_RCPP
}
// get_gsm_pointer
SEXP get_gsm_pointer(arma::mat const& X, arma::mat const& XD, arma::mat const& Z, arma::vec const& y, double const eps, double const kappa, std::string const& link, unsigned const n_threads, arma::vec const& offset_eta, arma::vec const& offset_etaD);
RcppExport SEXP _survTMB_get_gsm_pointer(SEXP XSEXP, SEXP XDSEXP, SEXP ZSEXP, SEXP ySEXP, SEXP epsSEXP, SEXP kappaSEXP, SEXP linkSEXP, SEXP n_threadsSEXP, SEXP offset_etaSEXP, SEXP offset_etaDSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
  Rcpp::traits::input_parameter< arma::mat const& >::type XD(XDSEXP);
  Rcpp::traits::input_parameter< arma::mat const& >::type Z(ZSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
  Rcpp::traits::input_parameter< double const >::type eps(epsSEXP);
  Rcpp::traits::input_parameter< double const >::type kappa(kappaSEXP);
  Rcpp::traits::input_parameter< std::string const& >::type link(linkSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type n_threads(n_threadsSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type offset_eta(offset_etaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type offset_etaD(offset_etaDSEXP);
  rcpp_result_gen = Rcpp::wrap(get_gsm_pointer(X, XD, Z, y, eps, kappa, link, n_threads, offset_eta, offset_etaD));
  return rcpp_result_gen;
  END_RCPP
}
// gsm_eval_ll
double gsm_eval_ll(SEXP ptr, arma::vec const& beta, arma::vec const& gamma);
RcppExport SEXP _survTMB_gsm_eval_ll(SEXP ptrSEXP, SEXP betaSEXP, SEXP gammaSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type beta(betaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type gamma(gammaSEXP);
  rcpp_result_gen = Rcpp::wrap(gsm_eval_ll(ptr, beta, gamma));
  return rcpp_result_gen;
  END_RCPP
}
// gsm_eval_grad
arma::vec gsm_eval_grad(SEXP ptr, arma::vec const& beta, arma::vec const& gamma);
RcppExport SEXP _survTMB_gsm_eval_grad(SEXP ptrSEXP, SEXP betaSEXP, SEXP gammaSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type beta(betaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type gamma(gammaSEXP);
  rcpp_result_gen = Rcpp::wrap(gsm_eval_grad(ptr, beta, gamma));
  return rcpp_result_gen;
  END_RCPP
}
// gsm_eval_hess
arma::mat gsm_eval_hess(SEXP ptr, arma::vec const& beta, arma::vec const& gamma);
RcppExport SEXP _survTMB_gsm_eval_hess(SEXP ptrSEXP, SEXP betaSEXP, SEXP gammaSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type beta(betaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type gamma(gammaSEXP);
  rcpp_result_gen = Rcpp::wrap(gsm_eval_hess(ptr, beta, gamma));
  return rcpp_result_gen;
  END_RCPP
}
// get_pedigree_funcs
SEXP get_pedigree_funcs(Rcpp::List data, Rcpp::List parameters);
RcppExport SEXP _survTMB_get_pedigree_funcs(SEXP dataSEXP, SEXP parametersSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
  Rcpp::traits::input_parameter< Rcpp::List >::type parameters(parametersSEXP);
  rcpp_result_gen = Rcpp::wrap(get_pedigree_funcs(data, parameters));
  return rcpp_result_gen;
  END_RCPP
}
// pedigree_funcs_eval_lb
double pedigree_funcs_eval_lb(SEXP p, SEXP par);
RcppExport SEXP _survTMB_pedigree_funcs_eval_lb(SEXP pSEXP, SEXP parSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  rcpp_result_gen = Rcpp::wrap(pedigree_funcs_eval_lb(p, par));
  return rcpp_result_gen;
  END_RCPP
}
void pedigree_funcs_eval_grad(SEXP p, SEXP par, Rcpp::NumericVector out);
RcppExport SEXP _survTMB_pedigree_funcs_eval_grad(SEXP pSEXP, SEXP parSEXP, SEXP outSEXP) {
  BEGIN_RCPP
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type out(outSEXP);
  pedigree_funcs_eval_grad(p, par, out);
  return R_NilValue;
  END_RCPP
}
// joint_start_ll
arma::vec joint_start_ll(arma::vec const& Y, arma::vec const& tstart, arma::vec const& tstop, arma::vec const& omega, arma::vec const& delta, arma::mat const& Z, unsigned const n_nodes, arma::vec const& coefs, bool const grad, bool const use_log, std::string const basis_type);
RcppExport SEXP _survTMB_joint_start_ll(SEXP YSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP omegaSEXP, SEXP deltaSEXP, SEXP ZSEXP, SEXP n_nodesSEXP, SEXP coefsSEXP, SEXP gradSEXP, SEXP use_logSEXP, SEXP basis_typeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< arma::vec const& >::type Y(YSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type tstart(tstartSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type tstop(tstopSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type omega(omegaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type delta(deltaSEXP);
  Rcpp::traits::input_parameter< arma::mat const& >::type Z(ZSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type n_nodes(n_nodesSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type coefs(coefsSEXP);
  Rcpp::traits::input_parameter< bool const >::type grad(gradSEXP);
  Rcpp::traits::input_parameter< bool const >::type use_log(use_logSEXP);
  Rcpp::traits::input_parameter< std::string const >::type basis_type(basis_typeSEXP);
  rcpp_result_gen = Rcpp::wrap(joint_start_ll(Y, tstart, tstop, omega, delta, Z, n_nodes, coefs, grad, use_log, basis_type));
  return rcpp_result_gen;
  END_RCPP
}
// get_joint_funcs
SEXP get_joint_funcs(Rcpp::List data, Rcpp::List parameters);
RcppExport SEXP _survTMB_get_joint_funcs(SEXP dataSEXP, SEXP parametersSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
  Rcpp::traits::input_parameter< Rcpp::List >::type parameters(parametersSEXP);
  rcpp_result_gen = Rcpp::wrap(get_joint_funcs(data, parameters));
  return rcpp_result_gen;
  END_RCPP
}
// joint_funcs_eval_lb
double joint_funcs_eval_lb(SEXP p, SEXP par);
RcppExport SEXP _survTMB_joint_funcs_eval_lb(SEXP pSEXP, SEXP parSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  rcpp_result_gen = Rcpp::wrap(joint_funcs_eval_lb(p, par));
  return rcpp_result_gen;
  END_RCPP
}
// joint_funcs_eval_grad
Rcpp::NumericVector joint_funcs_eval_grad(SEXP p, SEXP par);
RcppExport SEXP _survTMB_joint_funcs_eval_grad(SEXP pSEXP, SEXP parSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  Rcpp::traits::input_parameter< SEXP >::type par(parSEXP);
  rcpp_result_gen = Rcpp::wrap(joint_funcs_eval_grad(p, par));
  return rcpp_result_gen;
  END_RCPP
}
// get_orth_poly
List get_orth_poly(arma::vec const& x, unsigned const degree);
RcppExport SEXP _survTMB_get_orth_poly(SEXP xSEXP, SEXP degreeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< arma::vec const& >::type x(xSEXP);
  Rcpp::traits::input_parameter< unsigned const >::type degree(degreeSEXP);
  rcpp_result_gen = Rcpp::wrap(get_orth_poly(x, degree));
  return rcpp_result_gen;
  END_RCPP
}
// predict_orth_poly
arma::mat predict_orth_poly(arma::vec const& x, arma::vec const& alpha, arma::vec const& norm2);
RcppExport SEXP _survTMB_predict_orth_poly(SEXP xSEXP, SEXP alphaSEXP, SEXP norm2SEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< arma::vec const& >::type x(xSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type alpha(alphaSEXP);
  Rcpp::traits::input_parameter< arma::vec const& >::type norm2(norm2SEXP);
  rcpp_result_gen = Rcpp::wrap(predict_orth_poly(x, alpha, norm2));
  return rcpp_result_gen;
  END_RCPP
}
// fix_atomic_seqfault
void fix_atomic_seqfault();
RcppExport SEXP _survTMB_fix_atomic_seqfault() {
  BEGIN_RCPP
  fix_atomic_seqfault();
  return R_NilValue;
  END_RCPP
}
// setup_atomic_cache
void setup_atomic_cache(size_t const n_nodes, std::string const type, std::string const link);
RcppExport SEXP _survTMB_setup_atomic_cache(SEXP n_nodesSEXP, SEXP typeSEXP, SEXP linkSEXP) {
  BEGIN_RCPP
  Rcpp::traits::input_parameter< size_t const >::type n_nodes(n_nodesSEXP);
  Rcpp::traits::input_parameter< std::string const >::type type(typeSEXP);
  Rcpp::traits::input_parameter< std::string const >::type link(linkSEXP);
  setup_atomic_cache(n_nodes, type, link);
  return R_NilValue;
  END_RCPP
}
// pedigree_get_size
Rcpp::List pedigree_get_size(SEXP p);
RcppExport SEXP _survTMB_pedigree_get_size(SEXP pSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< SEXP >::type p(pSEXP);
  rcpp_result_gen = Rcpp::wrap(pedigree_get_size(p));
  return rcpp_result_gen;
  END_RCPP
}
// clear_cppad_mem
int clear_cppad_mem(unsigned const max_n_threads, bool const keep_work_space);
RcppExport SEXP _survTMB_clear_cppad_mem(SEXP max_n_threadsSEXP, SEXP keep_work_spaceSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< unsigned const >::type max_n_threads(max_n_threadsSEXP);
  Rcpp::traits::input_parameter< bool const >::type keep_work_space(keep_work_spaceSEXP);
  rcpp_result_gen = Rcpp::wrap(clear_cppad_mem(max_n_threads, keep_work_space));
  return rcpp_result_gen;
  END_RCPP
}
// set_n_threads
int set_n_threads(int const n_threads);
RcppExport SEXP _survTMB_set_n_threads(SEXP n_threadsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< int const >::type n_threads(n_threadsSEXP);
  rcpp_result_gen = Rcpp::wrap(set_n_threads(n_threads));
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
  {"_survTMB_joint_start_ll", (DL_FUNC) &_survTMB_joint_start_ll, 11},
  {"_survTMB_get_joint_funcs", (DL_FUNC) &_survTMB_get_joint_funcs, 2},
  {"_survTMB_joint_funcs_eval_lb", (DL_FUNC) &_survTMB_joint_funcs_eval_lb, 2},
  {"_survTMB_joint_funcs_eval_grad", (DL_FUNC) &_survTMB_joint_funcs_eval_grad, 2},
  {"_survTMB_get_commutation", (DL_FUNC) &_survTMB_get_commutation, 2},
  {"_survTMB_get_gsm_pointer", (DL_FUNC) &_survTMB_get_gsm_pointer, 10},
  {"_survTMB_gsm_eval_ll", (DL_FUNC) &_survTMB_gsm_eval_ll, 3},
  {"_survTMB_gsm_eval_grad", (DL_FUNC) &_survTMB_gsm_eval_grad, 3},
  {"_survTMB_gsm_eval_hess", (DL_FUNC) &_survTMB_gsm_eval_hess, 3},
  {"_survTMB_get_pedigree_funcs", (DL_FUNC) &_survTMB_get_pedigree_funcs, 2},
  {"_survTMB_pedigree_funcs_eval_lb", (DL_FUNC) &_survTMB_pedigree_funcs_eval_lb, 2},
  {"_survTMB_pedigree_funcs_eval_grad", (DL_FUNC) &_survTMB_pedigree_funcs_eval_grad, 3},
  {"_survTMB_get_orth_poly", (DL_FUNC) &_survTMB_get_orth_poly, 2},
  {"_survTMB_predict_orth_poly", (DL_FUNC) &_survTMB_predict_orth_poly, 3},
  {"_survTMB_fix_atomic_seqfault", (DL_FUNC) &_survTMB_fix_atomic_seqfault, 0},
  {"_survTMB_setup_atomic_cache", (DL_FUNC) &_survTMB_setup_atomic_cache, 3},
  {"_survTMB_pedigree_get_size", (DL_FUNC) &_survTMB_pedigree_get_size, 1},
  {"_survTMB_clear_cppad_mem", (DL_FUNC) &_survTMB_clear_cppad_mem, 2},
  {"_survTMB_set_n_threads", (DL_FUNC) &_survTMB_set_n_threads, 1},
  {"_survTMB_psqn_get_mgsm_funcs", (DL_FUNC) &_survTMB_psqn_get_mgsm_funcs, 9},
  {"_survTMB_psqn_optim_mgsm", (DL_FUNC) &_survTMB_psqn_optim_mgsm, 11},
  {"_survTMB_psqn_optim_mgsm_private", (DL_FUNC) &_survTMB_psqn_optim_mgsm_private, 7},
  {"_survTMB_eval_psqn_mgsm", (DL_FUNC) &_survTMB_eval_psqn_mgsm, 3},
  {"_survTMB_grad_psqn_mgsm", (DL_FUNC) &_survTMB_grad_psqn_mgsm, 3},
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
