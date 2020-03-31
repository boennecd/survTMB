#include "gva.h"
#include "gva-utils.h"
#include "utils.h"

using namespace survTMB;
using GaussHermite::GVA::mlogit_integral;
using GaussHermite::GVA::probit_integral;

namespace {

/* util class and functions to compute the conditional density of the
 * observed outcomes */

template<class Type>
struct GVA_cond_dens_data {
  /* input dependend objects */
  Type const eps,
         eps_log = Type(log(eps)),
           kappa;
  unsigned const n_nodes;

  /* potentially needed constants */
  Type const mlog_2_pi_half = Type(-log(2 * M_PI) / 2.),
                        two = Type(2.);

  GVA_cond_dens_data
    (Type const &eps, Type const &kappa, unsigned const n_nodes):
    eps(eps), kappa(kappa), n_nodes(n_nodes) { }
};

#define GVA_COND_DENS_ARGS                                     \
  Type const &eta_fix, Type const &etaD_fix,                   \
  Type const &event, Type const &va_mean, Type const &va_sd,   \
  Type const &va_var

/* computes the conditional density term for the PH (log-log) link
 * function */
template<class Type>
struct ph final : public GVA_cond_dens_data<Type> {
  using GVA_cond_dens_data<Type>::GVA_cond_dens_data;

  Type operator()(GVA_COND_DENS_ARGS) const {
    Type const eta = eta_fix  + va_mean,
                h = etaD_fix * exp(eta),
                H = exp(eta + va_var / this->two),
           if_low = event * this->eps_log - H - h * h * this->kappa,
           if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, this->eps, if_ok, if_low);
  }
};

/* computes the conditional density term for the PO (-logit) link
 * function */
template<class Type>
struct po final : public GVA_cond_dens_data<Type> {
  using GVA_cond_dens_data<Type>::GVA_cond_dens_data;

  Type operator()(GVA_COND_DENS_ARGS) const {
    Type const eta = eta_fix + va_mean,
                 H = mlogit_integral(
                   va_mean, va_sd, eta_fix, this->n_nodes),
                 h = etaD_fix * exp(eta - H),
            if_low = event * this->eps_log - H - h * h * this->kappa,
            if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, this->eps, if_ok, if_low);
  }
};

/* computes the conditional density term for the probit (-probit) link
 * function */
template<class Type>
struct probit : public GVA_cond_dens_data<Type> {
  using GVA_cond_dens_data<Type>::GVA_cond_dens_data;

  Type operator()(GVA_COND_DENS_ARGS) const {
    Type const H = probit_integral(va_mean, va_sd, -eta_fix, this->n_nodes),
            diff = eta_fix + va_mean,
               h = etaD_fix * exp(
                 this->mlog_2_pi_half - diff * diff / this->two -
                   va_var / this->two + H),
          if_low = event * this->eps_log - H - h * h * this->kappa,
          if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, this->eps, if_ok, if_low);
  }
};

#undef GVA_COND_DENS_ARGS

template<class Type>
void GVA_comp(COMMON_ARGS(Type), vector<Type> const &theta_VA,
              unsigned const n_nodes, unsigned const rng_dim){
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
  auto main_loop = [&](auto const &func){
    unsigned i = 0;
    for(unsigned g = 0; g < grp_size.size(); ++g){
      unsigned const n_members = grp_size[g];
      /* is this our cluster? */
      if(!is_my_region(*result.obj)){
        i += n_members;
        result.obj->parallel_region();
        continue;
      }

      /* get VA parameters */
      vecT const &va_mu          = va_means[g];
      matrix<Type> const &va_var = va_vcovs[g];

      /* compute conditional density terms from outcomes */
      unsigned const end = n_members + i;
      Type terms(0);
      for(; i < end; ++i){
        vecT const z = Z.row(i);
        Type const err_mean = (z * va_mu).sum(),
                   err_var  = (z * vecT(va_var * z)).sum(),
                   err_sd   = sqrt(err_var);

        terms += func(
          eta_fix[i], etaD_fix[i], event[i], err_mean, err_sd, err_var);
      }

      result -= terms;
    }
  };

  if(link == "PH"){
    ph<Type> const func(eps, kappa, n_nodes);
    main_loop(func);

  } else if (link == "PO"){
    po<Type> const func(eps, kappa, n_nodes);
    main_loop(func);

  } else if (link == "probit"){
    probit<Type> const func(eps, kappa, n_nodes);
    main_loop(func);

  } else
    error("'%s' not implemented", link.c_str());

  if(!is_my_region(*result.obj))
    /* only have to add one more term so just return */
    return;

  /* handle terms from random effect log density and VA log density */
  Type lb_term(0.);
  {
    matrix<Type> va_cov_sum(rng_dim, rng_dim);
    va_cov_sum.setZero();

    for(unsigned g = 0; g < n_groups; ++g){
      lb_term += atomic::logdet(va_vcovs[g]) -
        (va_means[g] * vecT(vcov_inv * va_means[g])).sum();
      va_cov_sum += va_vcovs[g];

    }
    matrix<Type> const mat_prod = va_cov_sum * vcov_inv;
    lb_term -= mat_prod.trace();
  }

  /* add final log determinant term and constants */
  lb_term += Type(n_groups) * (-log_det_vcov + Type(rng_dim));
  result -= lb_term / Type(2.);
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

  return GVA_comp(COMMON_CALL, theta_VA, n_nodes, rng_dim);
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

} // namespace survTMB
