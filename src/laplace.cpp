#include "laplace.h"
#include "dmvnorm_log.h"
#include "utils.h"
#include "pnorm_log.h"
#include "utils.h"

namespace survTMB {

template<class Type>
Type laplace_PH_terms
  (Type const &eta, Type const &etaD, Type const &event,
   Type const &eps, Type const &eps_log, Type const &kappa){
  Type const H = exp(eta),
             h = etaD * H,
        if_low = event * eps_log - H - h * h * kappa,
        if_ok  = event * log(h)  - H;
  return CppAD::CondExpGe(h, eps, if_ok, if_low);
}

template<class Type>
Type laplace_PO_terms
  (Type const &eta, Type const &etaD, Type const &event,
   Type const &eps, Type const &eps_log, Type const &kappa){
  Type const too_large(30.),
                   one(1.);

  Type const H = CppAD::CondExpGe(
    eta, too_large, eta, log(one + exp(eta))),
             h = etaD * exp(eta - H),
        if_low = event * eps_log - H - h * h * kappa,
        if_ok  = event * log(h)  - H;
  return CppAD::CondExpGe(h, eps, if_ok, if_low);
}

template<class Type>
Type laplace_probit_terms
  (Type const &eta, Type const &etaD, Type const &event,
   Type const &eps, Type const &eps_log, Type const &kappa){
  Type const tiny(std::numeric_limits<double>::epsilon()),
             zero(0.),
              one(1.);

  Type const H = -pnorm_log(-eta),
             h = etaD * dnorm(-eta, zero, one) /
               (pnorm(-eta) + tiny),
        if_low = event * eps_log - H - h * h * kappa,
        if_ok  = event * log(h)  - H;
  return CppAD::CondExpGe(h, eps, if_ok, if_low);
}

template<class Type>
void laplace(COMMON_ARGS(Type), matrix<Type> const &u){
  /* checks */
  unsigned const rng_dim = survTMB::get_rng_dim(theta);
  {
    int const *max_grp = std::max_element(grp.data(), grp.data() + grp.size());
    if(!max_grp or *max_grp != u.cols() - 1L)
      error("Invalid 'grp' or 'u'");
    if(rng_dim != u.rows())
      error("Invalid 'u'");
  }

  /* log-likelihood terms from conditional distribution of the observed
   * outcomes */
  Type  const eps_log = log(eps);
  auto const b_vec = b.matrix().col(0);

  /* compute terms from conditional density */
  typedef Type (*loop_func)(
      Type const&, Type const&, Type const&,
      Type const&, Type const&, Type const&);

  auto cond_dens_loop = [&](loop_func func){
    unsigned i(0L);
    for(unsigned g = 0; g < grp_size.size(); ++g){
      unsigned const n_members = grp_size[g];
      /* do we need to anything on this thread? */
      if(!is_my_region(*result.obj)){
        i += n_members;
        result.obj->parallel_region();
        continue;
      }

      /* compute linear predictor etc. */
      auto const u_g = u.col(g);
      unsigned const end_i(i + n_members);
      vector<Type> const eta = ([&](){
        vector<Type> out(n_members);
        unsigned k(0L);
        for(unsigned j = i; j < end_i; ++j, ++k){
          out[k]  = (X.row(j) * b_vec)[0];
          out[k] += Z.row(j) * u_g;
        }

        return out;
      })();
      vector<Type> const etaD = ([&](){
        vector<Type> out(n_members);
        unsigned k(0L);
        for(unsigned j = i; j < end_i; ++j, ++k)
          out[k] = (XD.row(j) * b_vec)[0];

        return out;
      })();

      Type next_term(0.);
      for(unsigned k = 0; k < n_members; ++k, ++i)
        next_term += func(
          eta[k], etaD[k], event[i], eps, eps_log, kappa);

      result -= next_term;
    }
  };

  if(link == "PH")
    cond_dens_loop(laplace_PH_terms<Type>);
  else if (link == "PO")
    cond_dens_loop(laplace_PO_terms<Type>);
  else if(link == "probit")
    cond_dens_loop(laplace_probit_terms<Type>);
  else
    error("'%s' not implemented", link.c_str());

  /* log-likelihood terms from random effect density */
  result -= mult_var_dens(theta, u);
}

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;
using namespace survTMB;

template void laplace<double>(
    COMMON_ARGS(double), matrix<double> const &u);
template void laplace<ADd>(
    COMMON_ARGS(ADd), matrix<ADd> const &u);
template void laplace<ADdd>(
    COMMON_ARGS(ADdd), matrix<ADdd> const &u);
template void laplace<ADddd>(
    COMMON_ARGS(ADddd), matrix<ADddd> const &u);

} // namespace survTMB


