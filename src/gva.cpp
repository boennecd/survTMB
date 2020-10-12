#include "gva.h"
#include "gva-utils.h"
#include "utils.h"

using namespace survTMB;
using namespace GaussHermite::GVA;

namespace {
template<class Type, template <class> class Accumlator>
void GVA_comp(COMMON_ARGS(Type, Accumlator), vector<Type> const &theta_VA,
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
  bool const is_in_parallel = CppAD::thread_alloc::in_parallel();
#define MAIN_LOOP(func)                                        \
  {                                                            \
    unsigned i = 0;                                            \
    for(unsigned g = 0; g < grp_size.size(); ++g){             \
      unsigned const n_members = grp_size[g];                  \
      /* is this our cluster? */                               \
      if(is_in_parallel and !is_my_region(*result.obj)){       \
        i += n_members;                                        \
        result.obj->parallel_region();                         \
        continue;                                              \
      }                                                        \
                                                               \
      /* get VA parameters */                                  \
      vecT const &va_mu          = va_means[g];                \
      matrix<Type> const &va_var = va_vcovs[g];                \
                                                               \
      /* compute conditional density terms from outcomes */    \
      unsigned const end = n_members + i;                      \
      Type terms(0);                                           \
      for(; i < end; ++i){                                     \
        vecT const z = Z.row(i);                               \
        Type const err_mean = vec_dot(z, va_mu),               \
                   err_var  = quad_form_sym(z, va_var),        \
                   err_sd   = sqrt(err_var);                   \
                                                               \
        terms += func(                                         \
          eta_fix[i], etaD_fix[i], event[i], err_mean, err_sd, \
          err_var);                                            \
      }                                                        \
                                                               \
      result -= terms;                                         \
    }                                                          \
  }

  if(link == "PH"){
    ph<Type> const func(eps, kappa, n_nodes);
    MAIN_LOOP(func);

  } else if (link == "PO"){
    po<Type> const func(eps, kappa, n_nodes);
    MAIN_LOOP(func);

  } else if (link == "probit"){
    probit<Type> const func(eps, kappa, n_nodes);
    MAIN_LOOP(func);

  } else
    error("'%s' not implemented", link.c_str());

#undef MAIN_LOOP

  if(!is_my_region(*result.obj))
    /* only have to add one more term so just return */
    return;

  /* handle terms from random effect log density and VA log density */
  Type lb_term(0.);
  {
    matrix<Type> va_cov_sum(rng_dim, rng_dim);
    va_cov_sum.setZero();

    for(unsigned g = 0; g < n_groups; ++g){
      lb_term +=
        atomic::logdet(va_vcovs[g]) - quad_form_sym(va_means[g], vcov_inv);
      va_cov_sum += va_vcovs[g];

    }

    lb_term -= mat_mult_trace(va_cov_sum, vcov_inv);
  }

  /* add final log determinant term and constants */
  lb_term += Type(n_groups) * (-log_det_vcov + Type(rng_dim));
  result -= lb_term / Type(2.);
}

} // namespace

namespace survTMB {

template<class Type, template <class> class Accumlator>
void GVA(COMMON_ARGS(Type, Accumlator), vector<Type> const &theta_VA,
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

#define DEF_GVA(TYPE, ACCUMLATOR)                              \
  template void GVA<TYPE, ACCUMLATOR>                          \
    (COMMON_ARGS(TYPE, ACCUMLATOR),                            \
     vector<TYPE> const &theta_VA,                             \
     unsigned const n_nodes)

DEF_GVA(double, parallel_accumulator);
DEF_GVA(ADd   , parallel_accumulator);
DEF_GVA(ADdd  , parallel_accumulator);
DEF_GVA(ADddd , parallel_accumulator);

DEF_GVA(double, accumulator_mock);
DEF_GVA(ADd   , accumulator_mock);
DEF_GVA(ADdd  , accumulator_mock);
DEF_GVA(ADddd , accumulator_mock);

} // namespace survTMB
