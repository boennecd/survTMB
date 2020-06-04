#include "snva.h"
#include "snva-utils.h"
#include <math.h>

namespace {

using namespace GaussHermite;
using namespace GaussHermite::SNVA;
using namespace survTMB;

template<class Type, template <class> class Accumlator>
void SNVA_comp
  (COMMON_ARGS(Type, Accumlator), vector<Type> const &theta_VA,
   unsigned const n_nodes, std::string const &param_type,
   unsigned const rng_dim){
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
  std::vector<vecT > va_mus,
                    va_rhos;
  std::vector<matrix<Type> > va_lambdas;

#define SET_PARAMS(meth_use)                                   \
  auto const input = meth_use(                                 \
    &theta_VA[0], theta_VA.size(), rng_dim);                   \
  va_mus     = move(input.va_mus);                             \
  va_rhos    = move(input.va_rhos);                            \
  va_lambdas = move(input.va_lambdas)

  if(param_type == "DP"){
    SET_PARAMS(SNVA_MD_theta_DP_to_DP);
  } else if(param_type == "CP_trans"){
    SET_PARAMS(SNVA_MD_theta_CP_trans_to_DP);
  } else
    error("SNVA_MD: param_type '%s' is not implemented",
          param_type.c_str());
#undef SET_PARAMS

  unsigned const n_groups = va_mus.size();

  /* assign constant and fixed effect objects */
  vecT const eta_fix = X  * b,
            etaD_fix = XD * b;
  Type const sqrt_2_pi(sqrt(M_2_PI)),
                   one(1.),
                   two(2.),
                 small(std::numeric_limits<double>::epsilon());

  /* assign object used in the variational distribution */
  std::vector<vecT> va_ds;
  va_ds.reserve(n_groups);
  for(unsigned g = 0; g < n_groups; ++g){
    vecT new_d = va_lambdas[g] * va_rhos[g];
    Type const denom = sqrt(one + vec_dot(va_rhos[g], new_d));
    new_d /= denom;
    va_ds.emplace_back(move(new_d));
  }

  /* handle terms from conditional density of observed outcomes */
  bool const is_in_parallel = CppAD::thread_alloc::in_parallel();
  ph    <Type>     ph_func(eps, kappa, n_nodes);
  po    <Type>     po_func(eps, kappa, n_nodes);
  probit<Type> probit_func(eps, kappa, n_nodes);
  {
    unsigned i = 0;
    for(unsigned g = 0; g < grp_size.size(); ++g){
      unsigned const n_members = grp_size[g];
      /* is this our cluster? */
      if(is_in_parallel and !is_my_region(*result.obj)){
        i += n_members;
        result.obj->parallel_region();
        continue;
      }

      vecT const &va_mu = va_mus[g],
                  &va_d = va_ds [g];
      matrix<Type> const &va_lambda = va_lambdas[g];

      Type term(0.);
      unsigned const end = n_members + i;
      for(; i < end; ++i){
        vecT const z = Z.row(i);

        Type const mu = vec_dot(z, va_mu),
                sd_sq = quad_form_sym(z, va_lambda),
                   sd = sqrt(sd_sq),
                    d = vec_dot(z, va_d),
                  rho = d / sd_sq / sqrt(one - d * d / sd_sq),
             d_scaled = sqrt_2_pi * d,
            dist_mean = mu + d_scaled,
             dist_var = sd_sq - d_scaled * d_scaled;

        if(link == "PH")
          term += ph_func(
            eta_fix[i], etaD_fix[i], event[i],
            mu, sd, rho, d, sd_sq, dist_mean, dist_var);
        else if(link == "PO")
          term += po_func(
            eta_fix[i], etaD_fix[i], event[i],
            mu, sd, rho, d, sd_sq, dist_mean, dist_var);
        else if(link == "probit")
          term += probit_func(
            eta_fix[i], etaD_fix[i], event[i],
            mu, sd, rho, d, sd_sq, dist_mean, dist_var);
        else
          error("'%s' not implemented", link.c_str());
      }

      result -= term;
    }
  }

  if(!is_my_region(*result.obj))
    /* only have to add one more term so just return */
    return;

  /* handle terms from random effect log density and VA log density */
  Type last_terms(0.);
  {
    matrix<Type> va_lambda_sum(rng_dim, rng_dim);
    va_lambda_sum.setZero();
    Type lb_t_mult_half(0.), lb_t_mult_other(0.);
    for(unsigned g = 0; g < n_groups; ++g){
      lb_t_mult_half +=
        atomic::logdet(va_lambdas[g]) - quad_form_sym(va_mus[g], vcov_inv);
      va_lambda_sum += va_lambdas[g];
      lb_t_mult_other -= quad_form(va_mus[g], vcov_inv, va_ds[g]);

      auto const llt_mat = va_lambdas[g].llt();
      vecT va_rho_scaled =
        (llt_mat.matrixU() * va_rhos[g].matrix()).array() + small;
      Type const r_L_r = vec_dot(va_rho_scaled, va_rho_scaled);
      last_terms -= entropy_term(r_L_r, n_nodes);

    }
    lb_t_mult_half -= mat_mult_trace(va_lambda_sum, vcov_inv);

    last_terms += lb_t_mult_half / two + lb_t_mult_other * sqrt_2_pi;
  }

  /* add final log determinant term and constants */
  last_terms += Type(n_groups) * (
    -log_det_vcov /  two + Type(rng_dim) /  two - Type(M_LN2));

  result -= last_terms;
}
} // namespace

namespace survTMB {

template<class Type, template <class> class Accumlator>
void SNVA(COMMON_ARGS(Type, Accumlator), vector<Type> const &theta_VA,
          unsigned const n_nodes, std::string const &param_type){
  /* checks */
  unsigned const rng_dim = get_rng_dim(theta);
  {
    int const *max_grp =
      std::max_element(grp.data(), grp.data() + grp.size());
    /* require mean and full covariance matrix */
    unsigned const expe_size = (*max_grp + 1L) *
    (rng_dim * 2L + (rng_dim * (rng_dim + 1L)) / 2L);
    if(theta_VA.size() != expe_size)
      error("theta_VA.size(): %i. Expected %i.",
            theta_VA.size(), expe_size);
  }

  return SNVA_comp(COMMON_CALL, theta_VA, n_nodes, param_type, rng_dim);
}

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;

#define DEF_SNVA(TYPE, ACCUMLATOR)                             \
  template void SNVA<TYPE, ACCUMLATOR>                         \
    (COMMON_ARGS(TYPE, ACCUMLATOR),                            \
     vector<TYPE> const &theta_VA,                             \
     unsigned const n_nodes, std::string const &param_type)

DEF_SNVA(double, parallel_accumulator);
DEF_SNVA(ADd   , parallel_accumulator);
DEF_SNVA(ADdd  , parallel_accumulator);
DEF_SNVA(ADddd , parallel_accumulator);

DEF_SNVA(double, accumulator_mock);
DEF_SNVA(ADd   , accumulator_mock);
DEF_SNVA(ADdd  , accumulator_mock);
DEF_SNVA(ADddd , accumulator_mock);

} // namespace survTMB
