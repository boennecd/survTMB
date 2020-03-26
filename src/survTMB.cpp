#include "dmvnorm_log.h"
#include "utils.h"
#include "pnorm_log.h"

#include <cmath>
#include <algorithm>
#include <future>
#include <limits>

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace survTMB;

/* Common arguments used in approximations.
 *
 * Args:
 *   tobs: observed times.
 *   event: zero/one variables for whether the event is observed.
 *   X: fixed effect design matrix.
 *   XD: derivative of X wrt time.
 *   Z: random effect design matrix.
 *   grp: integer vector with group identifier.
 *   eps: small tolerance for positivity constraint on hazards.
 *   kappa: one-sided L2 penalty for positivity constraint.
 *   b: fixed effect coefficients.
 *   theta: covariance matrix parameters. See get_vcov_from_trian.
 */
#define COMMON_ARGS                                            \
  parallel_accumulator<Type> &result,                          \
  vector<Type> const &tobs, vector<Type> const &event,         \
  matrix<Type> const &X, matrix<Type> const &XD,               \
  matrix<Type> const &Z, vector<int> const &grp,               \
  Type const &eps, Type const &kappa, vector<Type> const &b,   \
  vector<Type> const &theta, std::string const &link,          \
  vector<int> &grp_size
#define COMMON_CALL                                            \
  result, tobs, event, X, XD, Z, grp, eps, kappa, b, theta,    \
  link, grp_size

namespace {

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

/* Computes the log-likelihood for given random effects.
 *
 * Args:
 *   u: [random effect dim] x [# groups] matrix with random effects.
 */
template<class Type>
void laplace(COMMON_ARGS, matrix<Type> const &u){
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
} // namespace

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* common data objects and parameters */
  DATA_STRING(app_type); // TODO: remove app_type
  DATA_VECTOR(tobs);

  DATA_VECTOR(event);
  DATA_MATRIX(X);
  DATA_MATRIX(XD);
  DATA_MATRIX(Z);
  DATA_IVECTOR(grp);
  DATA_STRING(link);
  DATA_IVECTOR(grp_size);
  DATA_INTEGER(n_threads)

  /* They are marked as parameter such that the user can change them
   * later */
  PARAMETER(eps);
  PARAMETER(kappa);

  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(theta);

  /* checks */
  {
    unsigned const n = tobs.size();
    auto check_rows = [n](matrix<Type> const &x, char const *msg){
      if(x.rows() != n)
        error(msg);
    };
    check_rows(X , "invalid 'X'");
    check_rows(XD, "invalid 'XD'");
    check_rows(Z , "invalid 'Z'");
    if(n != event.size())
      error("invalid 'event'");
    if(n != grp.size())
      error("invalid 'grp'");

    if(b.size()  != X.cols())
      error("invalid 'b'");
    if(XD.cols() != X.cols())
      error("invalid 'XD' (# columns)");

    if(Z.cols() != get_rng_dim(theta))
      error("invalid 'Z' (# columns: %d %d)", Z.cols(), get_rng_dim(theta));

    std::size_t grp_size_sum = 0;
    for(int i = 0; i < grp_size.size(); ++i)
      grp_size_sum += grp_size[i];
    if(n != grp_size_sum)
      error("invalid 'grp_size'");

    int g_max = 0;
    for(int i = 0; i < grp.size(); ++i)
      g_max = std::max(grp[i], g_max);
    if(g_max + 1L != grp_size.size())
      error("invalid 'grp_size' (does not match with number of groups)");

    for(int i = 1L; i < grp.size(); ++i){
      if(grp[i - 1L] > grp[i])
        error("'grp' is not sorted");
      if(grp[i - 1L] - grp[i] < -1L)
        error("too big gap in 'grp'");
    }
  }

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  parallel_accumulator<Type> result(this);

  /* perform approximations using method descibed at
   *   https://github.com/kaskr/adcomp/issues/233#issuecomment-306032192
   *
   * each branch may contain further parameters or data objects */
  if(app_type == "Laplace"){
    PARAMETER_MATRIX(u);
    laplace(COMMON_CALL, u);
    return result;

  }

  error("approximation method '%s' is not implemented", app_type.c_str());
  return Type(0);
}
