#ifndef LAPLACE_UTIL_H
#define LAPLACE_UTIL_H

#include "pnorm-log.h"

namespace survTMB {

template<class Type>
Type laplace_PH_terms
  (Type const &eta, Type const &etaD, Type const &event,
   Type const &eps, Type const &eps_log, Type const &kappa){
  Type const H = exp(eta),
             h = etaD         * H,
            hp = (etaD - eps) * H,
        if_low = event * (eps_log + eta) - H - hp * hp * kappa,
        if_ok  = event * log(h)          - H;
  return CppAD::CondExpGe(etaD, eps, if_ok, if_low);
}

template<class Type>
Type laplace_PO_terms
  (Type const &eta, Type const &etaD, Type const &event,
   Type const &eps, Type const &eps_log, Type const &kappa){
  static Type const too_large(30.),
                          one(1.);

  Type const H = CppAD::CondExpGe(
    eta, too_large, eta, log(one + exp(eta))),
            he = exp(eta - H),
             h = etaD         * he,
            hp = (etaD - eps) * he,
        if_low = event * (eps_log + eta - H) - H - hp * hp * kappa,
        if_ok  = event * log(h)              - H;
  return CppAD::CondExpGe(etaD, eps, if_ok, if_low);
}

template<class Type>
Type laplace_probit_terms
  (Type const &eta, Type const &etaD, Type const &event,
   Type const &eps, Type const &eps_log, Type const &kappa){
  static Type const half(.5),
            log_2pi_half(std::log(2 * M_PI) * .5);

  Type const H = -pnorm_log(-eta),
        hf_log = -log_2pi_half - eta * eta * half + H,
            hf = exp(hf_log),
            hp = (etaD - eps) * hf,
        if_low = event * (eps_log   + hf_log) - H - hp * hp * kappa,
        if_ok  = event * (log(etaD) + hf_log) - H;
  return CppAD::CondExpGe(etaD, eps, if_ok, if_low);
}

} // namespace survTMB

#endif
