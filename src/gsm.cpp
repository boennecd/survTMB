#include "gsm.h"

namespace gsm_objs {
double gsm_probit::g_log() const {
  return pnrm_log;
}
double gsm_probit::gp() const{
  return -exp(dnrm_log);
}
double gsm_probit::gp_g() const{
  return -exp(dnrm_log - pnrm_log);
}
double gsm_probit::gpp_gp() const {
  return -eta;
}
double gsm_probit::gpp() const {
  return eta * exp(dnrm_log);
}
double gsm_probit::d_gp_g() const {
  double const log_gp_g = dnrm_log - pnrm_log,
                gpp_g = eta * exp(log_gp_g),
              gp_sq_g_sq = exp(2 * log_gp_g);
  return gpp_g - gp_sq_g_sq;
}
double gsm_probit::d_gpp_gp() const {
  return -1.;
}


double gsm_ph::g_log() const {
  return -exp_eta;
}
double gsm_ph::gp() const {
  return -exp(eta - exp_eta);
}
double gsm_ph::gp_g() const {
  return -exp_eta;
}
double gsm_ph::gpp_gp() const {
  return 1. - exp_eta;
}
double gsm_ph::gpp() const {
  return gp() * (1 - exp_eta);
}
double gsm_ph::d_gp_g() const {
  return -exp_eta;
}
double gsm_ph::d_gpp_gp() const {
  return -exp_eta;
}


constexpr double const logit_too_large = 30;

double gsm_logit::g_log() const {
  return eta > logit_too_large ? -eta : -log(exp_eta_p1);
}
double gsm_logit::gp() const {
  return eta > logit_too_large ? 0. :  -exp_eta / exp_eta_p1 / exp_eta_p1;
}
double gsm_logit::gp_g() const {
  return eta > logit_too_large ? -1. : -exp_eta / exp_eta_p1;
}
double gsm_logit::gpp_gp() const {
  return eta > logit_too_large ? -1. : -(exp_eta - 1.) / exp_eta_p1;
}
double gsm_logit::gpp() const {
  return eta > logit_too_large ?
    0. : exp_eta * (exp_eta - 1.) / exp_eta_p1 / exp_eta_p1 / exp_eta_p1;
}
double gsm_logit::d_gp_g() const {
  return  eta > logit_too_large ? 0 : -exp_eta / exp_eta_p1 / exp_eta_p1;
}
double gsm_logit::d_gpp_gp() const {
  return eta > logit_too_large ? 0 : -2. * exp_eta / exp_eta_p1 / exp_eta_p1;
}

template<class T>
Rcpp::XPtr<gsm_base> create_gsm_obj(
    arma::mat const &X, arma::mat const &XD, arma::mat const &Z,
    arma::vec const &y, double const eps, double const kappa,
    unsigned const n_threads, arma::vec const &offset_eta,
    arma::vec const &offset_etaD){
  return Rcpp::XPtr<gsm_base>(new T(X, XD, Z, y, eps, kappa, n_threads,
                                    offset_eta, offset_etaD));
}
} // namespace gsm_objs

using namespace gsm_objs;

/** returns an XPtr to the abstract base class. */
// [[Rcpp::export(rng = false)]]
SEXP get_gsm_pointer(
    arma::mat const &X, arma::mat const &XD, arma::mat const &Z,
    arma::vec const &y, double const eps, double const kappa,
    std::string const &link, unsigned const n_threads,
    arma::vec const &offset_eta, arma::vec const &offset_etaD){
  if(link == "probit")
    return create_gsm_obj<gsm<gsm_probit> >(
        X, XD, Z, y, eps, kappa, n_threads, offset_eta, offset_etaD);
  else if(link == "PH")
    return create_gsm_obj<gsm<gsm_ph    > >(
        X, XD, Z, y, eps, kappa, n_threads, offset_eta, offset_etaD);
  else if(link == "PO")
    return create_gsm_obj<gsm<gsm_logit > >(
          X, XD, Z, y, eps, kappa, n_threads, offset_eta, offset_etaD);

  throw std::invalid_argument("get_gsm_pointer: link not implemented");
  return SEXP();
}

/** evaluates the log-likelihood. */
// [[Rcpp::export(rng = false)]]
double gsm_eval_ll(SEXP ptr, arma::vec const &beta, arma::vec const &gamma){
  Rcpp::XPtr<gsm_base> obj(ptr);
  return obj->log_likelihood(beta, gamma);
}

/** evaluates the gradient of the log-likelihood. */
// [[Rcpp::export(rng = false)]]
arma::vec gsm_eval_grad
  (SEXP ptr, arma::vec const &beta, arma::vec const &gamma){
  Rcpp::XPtr<gsm_base> obj(ptr);
  return obj->grad(beta, gamma);
}

/** evaluates the Hessian of the log-likelihood. */
// [[Rcpp::export(rng = false)]]
arma::mat gsm_eval_hess
  (SEXP ptr, arma::vec const &beta, arma::vec const &gamma){
  Rcpp::XPtr<gsm_base> obj(ptr);
  return obj->hess(beta, gamma);
}
