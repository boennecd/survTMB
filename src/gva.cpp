#include "gva.h"
#include "gva-utils.h"
#include "utils.h"

using namespace survTMB;
using GaussHermite::GVA::mlogit_integral;
using GaussHermite::GVA::probit_integral;

namespace {
#define GVA_COND_DENS_ARGS                                     \
  Type const &eta_fix, Type const &etaD_fix,                   \
  Type const &event, Type const &va_mean, Type const &va_sd,   \
  Type const &va_var

/* util class to compute the conditional density of the observed outcomes */
template<class Type>
class GVA_cond_dens {
  /* input dependend objects */
  Type const eps,
         eps_log = Type(log(eps)),
           kappa;
  GaussHermite::HermiteData<Type> const &GH_xw;

  /* potentially needed constants */
  Type const mlog_2_pi_half = Type(-log(2 * M_PI) / 2.),
                        two = Type(2.);
public:
  GVA_cond_dens
  (Type const &eps, Type const &kappa,
   GaussHermite::HermiteData<Type> const &GH_xw):
  eps(eps), kappa(kappa), GH_xw(GH_xw) { }

  /* computes the conditional density term for the PH (log-log) link
   * function */
  Type ph(GVA_COND_DENS_ARGS) const {
    Type const eta = eta_fix  + va_mean,
                 h = etaD_fix * exp(eta),
                 H = exp(eta + va_var / two),
            if_low = event * eps_log - H - h * h * kappa,
            if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }

  /* computes the conditional density term for the PO (-logit) link
   * function */
  Type po(GVA_COND_DENS_ARGS) const {
    Type const eta = eta_fix + va_mean,
                 H = mlogit_integral(
                   va_mean, va_sd, eta_fix, GH_xw),
                 h = etaD_fix * exp(eta - H),
            if_low = event * eps_log - H - h * h * kappa,
            if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }

  /* computes the conditional density term for the probit (-probit) link
   * function */
  Type probit(GVA_COND_DENS_ARGS) const {
    Type const H = probit_integral(
      va_mean, va_sd, -eta_fix, GH_xw),
            diff = eta_fix + va_mean,
               h = etaD_fix * exp(
                 mlog_2_pi_half - diff * diff / two -
                   va_var / two + H),
          if_low = event * eps_log - H - h * h * kappa,
          if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }
};

#undef GVA_COND_DENS_ARGS

template<class Type>
void GVA_comp(COMMON_ARGS(Type), vector<Type> const &theta_VA,
              GaussHermite::HermiteData<Type> const &GH_xw,
              unsigned const rng_dim){
  using Eigen::Dynamic;
  using vecT = vector<Type>;
  using std::move;

  /* get objects related to model covariance matrix */
  matrix<Type> const vcov = get_vcov_from_trian(theta);

  /* maybe not the best idea to matrix multiply by the precision matrix
   * instead of using solve... */
  Type log_det_vcov;
  matrix<Type> vcov_inv;
  vcov_inv = atomic::matinvpd(vcov, log_det_vcov);

  /* get objects from VA distribution */
  unsigned const  dt = (rng_dim * (rng_dim + 1L)) / 2L,
            n_groups = theta_VA.size() / (rng_dim + dt);
  std::vector<vecT >         va_means;
  std::vector<matrix<Type> > va_vcovs;
  va_means.reserve(n_groups);
  va_vcovs.reserve(n_groups);
  {
    Type const *t = &theta_VA[0];
    for(unsigned g = 0; g < n_groups; ++g){
      /* insert new mean vector */
      vecT mean_vec(rng_dim);
      for(unsigned i = 0; i < rng_dim; ++i)
        mean_vec[i] = *t++;
      va_means.emplace_back(move(mean_vec));

      /* insert new covariance matrix */
      va_vcovs.emplace_back(get_vcov_from_trian(t, rng_dim));
      t += dt;
    }
  }

  /* assign constant and fixed effects objects */
  vector<Type> const eta_fix  = X  * b,
                     etaD_fix = XD * b;

  /* handle terms from conditional density of observed outcomes */
  GVA_cond_dens<Type> const cond_dens(eps, kappa, GH_xw);

#define ADD_COND_DENS(func)                                    \
  {                                                            \
    unsigned i = 0;                                            \
    for(unsigned g = 0; g < grp_size.size(); ++g){             \
      vecT const &va_mu = va_means[g];                         \
      matrix<Type> const &va_var = va_vcovs[g];                \
                                                               \
      unsigned const end = grp_size[g] + i;                    \
      for(; i < end; ++i){                                     \
        vecT const z = Z.row(i);                               \
        Type const err_mean = (z * va_mu).sum(),               \
                   err_var  = (z * vecT(va_var * z)).sum(),    \
                   err_sd   = sqrt(err_var);                   \
                                                               \
        result -= cond_dens.func(                              \
          eta_fix[i], etaD_fix[i], event[i], err_mean,         \
          err_sd, err_var);                                    \
      }                                                        \
    }                                                          \
  }

  if(link == "PH")
    ADD_COND_DENS(ph)
  else if (link == "PO")
    ADD_COND_DENS(po)
  else if (link == "probit")
    ADD_COND_DENS(probit)
  else
    error("'%s' not implemented", link.c_str());

#undef ADD_COND_DENS

  /* handle terms from random effect log density and VA log density */
  {
    matrix<Type> va_cov_sum(rng_dim, rng_dim);
    va_cov_sum.setZero();
    Type lb_term(0.);
    for(unsigned g = 0; g < n_groups; ++g){
      lb_term += atomic::logdet(va_vcovs[g]) -
        (va_means[g] * vecT(vcov_inv * va_means[g])).sum();
      va_cov_sum += va_vcovs[g];

    }
    matrix<Type> const mat_prod = va_cov_sum * vcov_inv;
    lb_term -= mat_prod.trace();

    result -= lb_term / Type(2.);
  }

  /* add final log determinant term and constants */
  result -= Type(n_groups) * (-log_det_vcov + Type(rng_dim)) / Type(2.);
}

} // namespace

namespace survTMB {

template<class Type>
void GVA(COMMON_ARGS(Type), vector<Type> const &theta_VA,
         unsigned const n_nodes){
  unsigned const rng_dim = get_rng_dim(theta);
  {
    int const *max_grp =
      std::max_element(grp.data(), grp.data() + grp.size());
    /* require mean and full covariance matrix */
    unsigned const expe_size = (*max_grp + 1L) *
    (rng_dim + (rng_dim * (rng_dim + 1L)) / 2L);
    if(theta_VA.size() != expe_size)
      error("theta_VA.size(): %i. Expected %i.",
            theta_VA.size(), expe_size);
  }

  GaussHermite::HermiteData<Type> const
    GH_xw(GaussHermite::GaussHermiteData(n_nodes));
  vector<Type> GH_x(n_nodes), GH_w(n_nodes);

  return GVA_comp(COMMON_CALL, theta_VA, GH_xw, rng_dim);
}

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;

template void GVA<double>
  (COMMON_ARGS(double), vector<double> const &theta_VA,
   unsigned const n_nodes);
template void GVA<ADd>
  (COMMON_ARGS(ADd), vector<ADd> const &theta_VA,
   unsigned const n_nodes);
template void GVA<ADdd>
  (COMMON_ARGS(ADdd), vector<ADdd> const &theta_VA,
   unsigned const n_nodes);
template void GVA<ADddd>
  (COMMON_ARGS(ADddd), vector<ADddd> const &theta_VA,
   unsigned const n_nodes);

}
