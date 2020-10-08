#define INCLUDE_RCPP
#include "snva-utils.h"
#include "get-x.h"
#include "psqn.h"
#include "psqn-reporter.h"

using namespace GaussHermite::SNVA;
using size_use = int;

namespace {
class snva_psqn_func {
  mutable ADFun<double> this_ad_func;
  size_use const g_dim,
                 p_dim;

public:
  snva_psqn_func(Rcpp::List data, double const eps, double const kappa,
                 arma::vec const &b, arma::vec const &theta,
                 arma::vec const &theta_va, size_use const n_nodes,
                 std::string const &link):
  this_ad_func(([&](){
    // assign value needed to record the terms from the lower bound for the
    // cluster
    using Type = AD<double>;
    using vecAD = vector<Type>;

    DATA_VECTOR(tobs);
    DATA_VECTOR(event);
    DATA_MATRIX(X);
    DATA_MATRIX(XD);
    DATA_MATRIX(Z);

    // checks
    size_use const n = tobs.size();
    auto check_rows = [n](matrix<Type> const &x,
                          char const *msg){
      if(x.rows() != static_cast<int>(n))
        error(msg);
    };
    check_rows(X , "invalid 'X'");
    check_rows(XD, "invalid 'XD'");
    check_rows(Z , "invalid 'Z'");
    if(static_cast<int>(n) != event.size())
      error("invalid 'event'");

    size_use const p = b.n_elem;
    if(static_cast<int>(p)  != X.cols())
      error("invalid 'b'");
    if(static_cast<int>(p) != X.cols())
      error("invalid 'XD' (# columns)");

    size_use const rng_dim = survTMB::get_rng_dim(theta.n_elem);
    if(Z.cols() != static_cast<int>(rng_dim))
      error("invalid 'Z' (# columns: %d %d)", Z.cols(),
            survTMB::get_rng_dim(theta.n_elem));

    // record ADFun. First create the concatenated vector of model
    // parameters. Then do the copy back after having started the recording.
    ADFun<double> func_out;
    size_use const n_par = b.n_elem + theta.n_elem + theta_va.n_elem;
    vecAD par(n_par);
    {
      Type *pi = &par[0];
      auto add_ele = [&](arma::vec const &x){
        double const *xi = x.begin();
        for(size_use i = 0; i < static_cast<size_use>(x.n_elem);
            ++i, ++xi, ++pi)
          *pi = Type(*xi);
      };
      add_ele(b);
      add_ele(theta);
      add_ele(theta_va);
    }

    // start recording
    CppAD::Independent(par);

    vecAD b_ad       (b       .n_elem),
          theta_ad   (theta   .n_elem),
          theta_va_ad(theta_va.n_elem);
    {
      Type *pi = &par[0];
      auto fill_ad_vec = [&](vector<Type> &x){
        Type *xi = &x[0];
        for(int i = 0; i < x.size(); ++i, ++pi, ++xi)
          *xi = *pi;
      };
      fill_ad_vec(b_ad);
      fill_ad_vec(theta_ad);
      fill_ad_vec(theta_va_ad);
    }

    // do the computation
    auto va_par = SNVA_MD_theta_DP_to_DP(
      &theta_va_ad[0], theta_va_ad.size(), rng_dim);
    vecAD const va_mu  = std::move(va_par.va_mus[0]),
                va_rho = std::move(va_par.va_rhos[0]);
    matrix<Type> const va_lambda = std::move(va_par.va_lambdas[0]);

    /* assign constant and fixed effect objects */
    vecAD const eta_fix  = X  * b_ad,
                etaD_fix = XD * b_ad;
    Type const sqrt_2_pi(sqrt(M_2_PI)),
                    one(1.),
                    two(2.),
                  small(std::numeric_limits<double>::epsilon());

    /* assign object used in the variational distribution */
    vecAD va_d = va_lambda * va_rho;
    {
      Type const denom = sqrt(one + vec_dot(va_rho, va_d));
      va_d /= denom;
    }

    vector<Type> y(1);
    Type &term = y[0];
    term = Type(0.);
    /* handle terms from conditional density of observed outcomes */
    {
      ph    <Type>     ph_func(eps, kappa, n_nodes);
      po    <Type>     po_func(eps, kappa, n_nodes);
      probit<Type> probit_func(eps, kappa, n_nodes);

      for(size_use i = 0; i < n; ++i){
        vecAD const z = Z.row(i);

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
          throw std::invalid_argument(link + " not implemented");
      }
    }

    /* handle terms from random effect log density and VA log density */
    // TODO: the inversion is done for each term. This cannot be avoided
    //       as it is right now with CppAD but it does yield some overhead
    matrix<Type> const vcov = survTMB::get_vcov_from_trian(theta_ad);
    Type log_det_vcov;
    matrix<Type> vcov_inv;
    vcov_inv = atomic::matinvpd(vcov, log_det_vcov);

    Type half_term = atomic::logdet(va_lambda) -
      quad_form_sym(va_mu, vcov_inv) -
      mat_mult_trace(va_lambda, vcov_inv);
    half_term /= two;

    Type const entrop_arg = quad_form_sym(va_rho, va_lambda) + small;
    term += half_term - entropy_term(entrop_arg, n_nodes) -
      quad_form(va_mu, vcov_inv, va_d) * sqrt_2_pi;

    // add final log determinant term and constants
    term += -log_det_vcov / two + Type(rng_dim) / two - Type(M_LN2);

    // we work with the negative lower bound
    term *= -one;

    // stop recording and return
    func_out.Dependent(par, y);
    func_out.optimize();
    return func_out;
  })()),
  g_dim(b.n_elem + theta.n_elem),
  p_dim(theta_va.n_elem) { }

  snva_psqn_func(snva_psqn_func&& other):
    this_ad_func(std::move(other.this_ad_func)),
    g_dim(other.g_dim), p_dim(other.p_dim) { }

  snva_psqn_func(snva_psqn_func const &other):
    g_dim(other.g_dim), p_dim(other.p_dim) {
    // we may not use the copy constructor but we may use the assignment
    // constructor
    this_ad_func = other.this_ad_func;
  }

  size_t global_dim() const {
    return static_cast<size_t>(g_dim);
  }
  size_t private_dim() const {
    return static_cast<size_t>(p_dim);
  }

  double func(double const *point) const {
    // copy
    // TODO: avoid repeated memory allocation?
    vector<double> par(g_dim + p_dim);
    double *p = &par[0];
    for(size_use i = 0; i < g_dim + p_dim; ++i, ++p, ++point)
      *p = *point;

    // evaluate
    return this_ad_func.Forward(0, par)[0];
  }
  double grad
    (double const * __restrict__ point, double * __restrict__ gr) const {
    // copy
    // TODO: avoid repeated memory allocation?
    vector<double> par(g_dim + p_dim);
    double *p = &par[0];
    for(size_use i = 0; i < g_dim + p_dim; ++i, ++p, ++point)
      *p = *point;

    // evaluate
    double const out = this_ad_func.Forward(0, par)[0];
    vector<double> w(1);
    w[0] = 1;

    vector<double> grad = this_ad_func.Reverse(1, w);
    for(size_use i = 0; i < grad.size(); ++i, ++gr)
      *gr = grad[i];

    return out;
  }

  bool thread_safe() const {
    return true;
  }
};
} // namespace

using snva_psqn_optim = PSQN::optimizer<snva_psqn_func, PSQN::R_reporter>;

// [[Rcpp::export(rng = false)]]
SEXP psqn_get_mgsm_funcs
  (Rcpp::List data, double const eps, double const kappa,
   arma::vec const &b, arma::vec const &theta,
   arma::vec const &theta_va, int const n_nodes,
   std::string const &link, unsigned const max_threads){
  shut_up();
  setup_parallel_ad setup_ADd(max_threads);

  std::vector<snva_psqn_func> funcs;
  funcs.reserve(data.size());

  size_use const rng_dim = survTMB::get_rng_dim(theta.n_elem),
                n_grp_va = 2 * rng_dim + (rng_dim * (rng_dim + 1)) / 2;
  arma::vec theta_va_term(n_grp_va);
  size_use i(0);

  for(auto d : data){
    for(size_use j = 0; j < n_grp_va; ++j, ++i)
      theta_va_term[j] = theta_va[i];
    funcs.emplace_back(Rcpp::List(d), eps, kappa, b, theta, theta_va_term,
                       n_nodes, link);
  }

  return Rcpp::XPtr<snva_psqn_optim>
    (new snva_psqn_optim(funcs, max_threads));
}

// [[Rcpp::export(rng = false)]]
Rcpp::List psqn_optim_mgsm
  (Rcpp::NumericVector val, SEXP ptr, double const rel_eps,
   unsigned const max_it, unsigned const n_threads, double const c1 = 0.0001,
   double const c2 = .9, bool const use_bfgs = true, int const trace = 0L,
   double const cg_tol = .1, bool const strong_wolfe = true){
  Rcpp::XPtr<snva_psqn_optim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("psqn_optim_mgsm: invalid parameter size");

  Rcpp::NumericVector par = clone(val);
  optim->set_n_threads(n_threads);
  auto res = optim->optim(&par[0], rel_eps, max_it, c1, c2,
                          use_bfgs, trace, cg_tol, strong_wolfe);
  Rcpp::NumericVector counts = Rcpp::NumericVector::create(
    res.n_eval, res.n_grad,  res.n_cg);
  counts.names() = Rcpp::CharacterVector::create
    ("function", "gradient", "n_cg");

  return Rcpp::List::create(
    Rcpp::_["par"] = par, Rcpp::_["value"] = res.value,
    Rcpp::_["info"] = static_cast<int>(res.info),
    Rcpp::_["counts"] = counts,
    Rcpp::_["convergence"] = res.info == PSQN::info_code::converged);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector psqn_optim_mgsm_private
  (Rcpp::NumericVector val, SEXP ptr, double const rel_eps,
   unsigned const max_it, unsigned const n_threads, double const c1 = 0.0001,
   double const c2 = .9){
  Rcpp::XPtr<snva_psqn_optim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("psqn_optim_mgsm_private: invalid parameter size");

  Rcpp::NumericVector par = clone(val);
  optim->set_n_threads(n_threads);
  double const res = optim->optim_priv(&par[0], rel_eps, max_it, c1, c2);
  par.attr("value") = res;
  return par;
}

// [[Rcpp::export(rng = false)]]
double eval_psqn_mgsm(Rcpp::NumericVector val, SEXP ptr,
                      unsigned const n_threads){
  Rcpp::XPtr<snva_psqn_optim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("eval_psqn_mgsm: invalid parameter size");

  optim->set_n_threads(n_threads);
  return optim->eval(&val[0], nullptr, false);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector grad_psqn_mgsm(Rcpp::NumericVector val, SEXP ptr,
                                   unsigned const n_threads){
  Rcpp::XPtr<snva_psqn_optim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("grad_psqn_mgsm: invalid parameter size");

  Rcpp::NumericVector grad(val.size());
  optim->set_n_threads(n_threads);
  grad.attr("value") = optim->eval(&val[0], &grad[0], true);

  return grad;
}
