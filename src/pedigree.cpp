#define INCLUDE_RCPP
#include "get-x.h"
#include "snva-utils.h"
#include "gva-utils.h"
#include "psqn.h"
#include "psqn-reporter.h"
#include <vector>

namespace {
using namespace GaussHermite;

template<class Tout>
vector<Tout> get_args(arma::vec const &omega, arma::vec const &beta,
                      arma::vec const &log_sds, arma::vec const &va_par){
  vector<Tout> out(
      omega.n_elem + beta.n_elem + log_sds.n_elem + va_par.n_elem);
  Tout *o = &out[0];

  auto add_to_vec = [&](arma::vec const &x){
    for(size_t i = 0; i < x.n_elem; ++i, ++o)
      *o = Tout(x[i]);
  };
  add_to_vec(omega);
  add_to_vec(beta);
  add_to_vec(log_sds);
  add_to_vec(va_par);

  return out;
}

/* object to hold the data for a given cluster */
template<class Type>
struct cluster_data {
  Rcpp::List data;

  const DATA_VECTOR(event);
  const DATA_MATRIX(X);
  const DATA_MATRIX(XD);
  const DATA_MATRIX(Z);

  std::vector<matrix<Type> > const cor_mats = ([&](){
    Rcpp::List list_w_cor_mats = data["cor_mats"];

    std::vector<matrix<Type> > out;
    out.reserve(list_w_cor_mats.size());
    for(auto x : list_w_cor_mats)
      out.emplace_back(get_mat<Type>(x));

    return out;
  })();

  size_t const n_members = event.size();

  /** not thread-safe because of the R interaction! */
  cluster_data(Rcpp::List data): data(data) {
    if(static_cast<size_t>(event.size()) != n_members)
      throw std::invalid_argument("cluster_data<Type>: invalid event");
    else if(static_cast<size_t>(X.cols())  != n_members)
      throw std::invalid_argument("cluster_data<Type>: invalid X");
    else if(static_cast<size_t>(XD.cols()) != n_members)
      throw std::invalid_argument("cluster_data<Type>: invalid XD");
    else if(static_cast<size_t>(Z.cols()) != n_members)
      throw std::invalid_argument("cluster_data<Type>: invalid Z");
    for(auto &V : cor_mats)
      if(static_cast<size_t>(V.rows()) != n_members or
           static_cast<size_t>(V.cols()) != n_members)
        throw std::invalid_argument("cluster_data<Type>: invalid cor_mats");

    data = R_NilValue;
  }

  void check(int const n_nodes, size_t const n_priv,
             arma::vec const &omega, arma::vec const &beta,
             arma::vec const &log_sds, arma::vec const &va_par){
    if(n_nodes < 1L)
      throw std::invalid_argument("pedigree_element_func<Type>: invalid n_nodes");
    else if(static_cast<size_t>(va_par.size()) != n_priv)
      throw std::invalid_argument("pedigree_element_func<Type>: invalid va_par");
    if(static_cast<size_t>(X.rows()) != omega.size())
      throw std::invalid_argument("pedigree_element_func<Type>: invalid c_data (X)");
    else if(Z.rows() != beta.size())
      throw std::invalid_argument("pedigree_element_func<Type>: invalid c_data (Z)");
    else if(static_cast<size_t>(cor_mats.size()) != log_sds.size())
      throw std::invalid_argument(
          "pedigree_element_func<Type>: invalid c_data (cor_mats)");
    else if(static_cast<size_t>(XD.rows()) != omega.size())
      throw std::invalid_argument("pedigree_element_func<Type>: invalid c_data (XD)");
  }
};

class pedigree_element_func_gva  {
  using Type = AD<double>;
  using vecAD = vector<Type>;

  size_t const n_global,
               n_members,
               n_priv;
  mutable ADFun<double> this_ad_func;

public:
  /** Notice that va_par is specific to the element function. TODO: maybe
   *  not a good idea? */
  pedigree_element_func_gva
  (Rcpp::List data, int const n_nodes, std::string const &link,
   arma::vec const &omega, arma::vec const &beta, arma::vec const &log_sds,
   arma::vec const &va_par, double const eps_in, double const kappa_in):
  n_global(omega.size() + beta.size() + log_sds.size()),
  n_members(Rcpp::NumericVector(static_cast<SEXP>(data["event"])).size()),
  n_priv(n_members + (n_members * (n_members + 1L)) / 2L),
  this_ad_func(([&](){
    // setup data objects
    cluster_data<Type> c_dat(data);
    Type const eps(eps_in),
             kappa(kappa_in);

    auto &event    = c_dat.event;
    auto &X        = c_dat.X;
    auto &XD       = c_dat.XD;
    auto &Z        = c_dat.Z;
    auto &cor_mats = c_dat.cor_mats;

    c_dat.check(n_nodes, n_priv, omega, beta, log_sds, va_par);

    // record ADFun. First create the concatenated vector of model
    // parameters. Then do the copy back after having started the recording.
    vecAD par = get_args<Type>(omega, beta, log_sds, va_par);

    // start recording
    CppAD::Independent(par);

    // copy back
    vecAD omega_ad(omega  .n_elem),
          beta_ad (beta   .n_elem),
        log_sds_ad(log_sds.n_elem),
         va_par_ad(va_par .n_elem);
    {
      Type *pi = &par[0];
      auto fill_ad_vec = [&](vector<Type> &x){
        Type *xi = &x[0];
        for(int i = 0; i < x.size(); ++i, ++pi, ++xi)
          *xi = *pi;
      };
      fill_ad_vec(omega_ad);
      fill_ad_vec(beta_ad);
      fill_ad_vec(log_sds_ad);
      fill_ad_vec(va_par_ad);
    }

    // constants
    Type const one(1.),
               two(2.);

    vecAD const var_ad = ([&](){
      vecAD out(log_sds_ad.size());
      for(int i = 0; i < out.size(); ++i)
        out[i] = exp(two * log_sds_ad[i]);
      return out;
    })();

    // do the computation
    GVA::ph    <Type>     ph_func(eps, kappa, n_nodes);
    GVA::po    <Type>     po_func(eps, kappa, n_nodes);
    GVA::probit<Type> probit_func(eps, kappa, n_nodes);

    vecAD va_mu(n_members);
    for(size_t i = 0; i < n_members; ++i)
      va_mu[i] = va_par_ad[i];
    matrix<Type> const va_var = survTMB::get_vcov_from_trian
      (&va_par_ad[n_members], n_members);
    Type term(0);

    // handle terms from the conditional density
    for(size_t i = 0; i < n_members; ++i){
      vecAD const x  = X .col(i),
                  xd = XD.col(i),
                  z  = Z .col(i);

      Type const &err_var = va_var(i, i),
                   err_sd = sqrt(err_var),
                 err_mean = va_mu[i],
                  eta_fix = vec_dot(omega_ad, x) + vec_dot(beta_ad, z),
                 etaD_fix = vec_dot(omega_ad, xd);

      if(link == "PH")
        term += ph_func(
          eta_fix, etaD_fix, event[i], err_mean, err_sd, err_var);
      else if(link == "PO")
        term += po_func(
          eta_fix, etaD_fix, event[i], err_mean, err_sd, err_var);
      else if(link == "probit")
        term += probit_func(
          eta_fix, etaD_fix, event[i], err_mean, err_sd, err_var);
      else
        error("'%s' not implemented", link.c_str());
    }

    /* add prior and entropy terms */
    matrix<Type> sigma(n_members, n_members);
    sigma.setZero();
    for(int i = 0; i < var_ad.size(); ++i)
      for(size_t j = 0; j < n_members; ++j)
        for(size_t k = 0; k < n_members; ++k)
          sigma(k, j) += var_ad[i] * cor_mats[i](k, j);

    Type log_det_sigma;
    matrix<Type> const sigma_inv = atomic::matinvpd(sigma, log_det_sigma);

    term += (
      atomic::logdet(va_var) - quad_form_sym(va_mu, sigma_inv)
      - mat_mult_trace(va_var, sigma_inv) - log_det_sigma
      + Type(n_members)) / two;

    vector<Type> result(1L);
    result[0] = -term;

    // stop recording and return
    ADFun<double> func_out;
    func_out.Dependent(par, result);
    func_out.optimize();
    return func_out;
  })()) { }

  pedigree_element_func_gva(pedigree_element_func_gva&& other):
  n_global(other.n_global), n_members(other.n_members),
  n_priv(other.n_priv),
  this_ad_func(std::move(other.this_ad_func)){ }

  pedigree_element_func_gva(pedigree_element_func_gva const &other):
  n_global(other.n_global), n_members(other.n_members),
  n_priv(other.n_priv) {
    // we may not use the copy constructor but we may use the assignment
    // constructor
    this_ad_func = other.this_ad_func;
  }

  size_t global_dim() const {
    return n_global;
  }
  size_t private_dim() const {
    return n_priv;
  }

  double func(double const *point) const {
    // copy
    // TODO: avoid repeated memory allocation?
    vector<double> par(n_global + n_priv);
    double *p = &par[0];
    for(size_t i = 0; i < n_global + n_priv; ++i, ++p, ++point)
      *p = *point;

    // evaluate
    return this_ad_func.Forward(0, par)[0];
  }
  double grad
  (double const * __restrict__ point, double * __restrict__ gr) const {
    // copy
    // TODO: avoid repeated memory allocation?
    vector<double> par(n_global + n_priv);
    double *p = &par[0];
    for(size_t i = 0; i < n_global + n_priv; ++i, ++p, ++point)
      *p = *point;

    // evaluate
    double const out = this_ad_func.Forward(0, par)[0];
    vector<double> w(1);
    w[0] = 1;

    vector<double> grad = this_ad_func.Reverse(1, w);
    for(size_t i = 0; i < n_global + n_priv; ++i, ++gr)
      *gr = grad[i];

    return out;
  }

  bool thread_safe() const {
    return true;
  }
};

class pedigree_element_func_snva  {
  using Type = AD<double>;
  using vecAD = vector<Type>;

  size_t const n_global,
               n_members,
               n_priv;
  mutable ADFun<double> this_ad_func;

public:
  /** Notice that va_par is specific to the element function. TODO: maybe
   *  not a good idea? */
  pedigree_element_func_snva
  (Rcpp::List data, int const n_nodes, std::string const &link,
   arma::vec const &omega, arma::vec const &beta, arma::vec const &log_sds,
   arma::vec const &va_par, double const eps_in, double const kappa_in):
  n_global(omega.size() + beta.size() + log_sds.size()),
  n_members(Rcpp::NumericVector(static_cast<SEXP>(data["event"])).size()),
  n_priv(2L * n_members + (n_members * (n_members + 1L)) / 2L),
  this_ad_func(([&](){
    // setup data objects
    cluster_data<Type> c_dat(data);
    Type const eps(eps_in),
             kappa(kappa_in);

    auto &event    = c_dat.event;
    auto &X        = c_dat.X;
    auto &XD       = c_dat.XD;
    auto &Z        = c_dat.Z;
    auto &cor_mats = c_dat.cor_mats;

    c_dat.check(n_nodes, n_priv, omega, beta, log_sds, va_par);

    // record ADFun. First create the concatenated vector of model
    // parameters. Then do the copy back after having started the recording.
    vecAD par = get_args<Type>(omega, beta, log_sds, va_par);

    // start recording
    CppAD::Independent(par);

    // copy back
    vecAD omega_ad(omega  .n_elem),
          beta_ad (beta   .n_elem),
        log_sds_ad(log_sds.n_elem),
         va_par_ad(va_par .n_elem);
    {
      Type *pi = &par[0];
      auto fill_ad_vec = [&](vector<Type> &x){
        Type *xi = &x[0];
        for(int i = 0; i < x.size(); ++i, ++pi, ++xi)
          *xi = *pi;
      };
      fill_ad_vec(omega_ad);
      fill_ad_vec(beta_ad);
      fill_ad_vec(log_sds_ad);
      fill_ad_vec(va_par_ad);
    }

    // constants
    Type const sqrt_2_pi(sqrt(M_2_PI)),
                     one(1.),
                     two(2.),
              type_M_LN2(M_LN2);

    vector<Type> const var_ad = ([&](){
      vector<Type> out(log_sds_ad.size());
      for(int i = 0; i < out.size(); ++i)
        out[i] = exp(two * log_sds_ad[i]);
      return out;
    })();

    SNVA::SNVA_MD_input<Type> va_par_trans =
      SNVA::SNVA_MD_theta_DP_to_DP(
        va_par_ad.data(), va_par_ad.size(), n_members);
    matrix<Type> const &lambda = va_par_trans.va_lambdas[0];
    vecAD const &mu = va_par_trans.va_mus[0],
               &rho = va_par_trans.va_rhos[0];
    matrix<Type> const rho_mat = asMatrix(rho, rho.size(), 1L);

    // do the computation
    SNVA::ph    <Type>     ph_func(eps, kappa, n_nodes);
    SNVA::po    <Type>     po_func(eps, kappa, n_nodes);
    SNVA::probit<Type> probit_func(eps, kappa, n_nodes);

    vecAD const delta = ([&](){
      vector<Type> out = atomic::matmul(lambda, rho_mat);
      Type const denom = sqrt(one + vec_dot(out, rho));
      out /= denom;
      return out;
    })();
    Type term(0.);

    // handle terms from the conditional density
    for(size_t i = 0; i < n_members; ++i){
      vecAD const x  = X .col(i),
                  xd = XD.col(i),
                  z  = Z .col(i);

      Type const &sd_sq = lambda(i, i),
                     sd = sqrt(sd_sq),
                     &d = delta[i],
                  rho_i = d / sd_sq / sqrt(one - d * d / sd_sq),
               d_scaled = sqrt_2_pi * d,
              dist_mean = mu[i] + d_scaled,
               dist_var = sd_sq - d_scaled * d_scaled,
                eta_fix = vec_dot(omega_ad, x) + vec_dot(beta_ad, z),
               etaD_fix = vec_dot(omega_ad, xd);

      if(link == "PH")
        term += ph_func(
          eta_fix, etaD_fix, event[i], mu[i], sd, rho_i, d, sd_sq,
          dist_mean, dist_var);
      else if(link == "PO")
        term += po_func(
          eta_fix, etaD_fix, event[i], mu[i], sd, rho_i, d, sd_sq,
          dist_mean, dist_var);
      else if(link == "probit")
        term += probit_func(
          eta_fix, etaD_fix, event[i], mu[i], sd, rho_i, d, sd_sq,
          dist_mean, dist_var);
      else
        error("'%s' not implemented", link.c_str());
    }

    /* add prior and entropy terms */
    matrix<Type> sigma(n_members, n_members);
    sigma.setZero();
    for(int i = 0; i < var_ad.size(); ++i)
      for(size_t j = 0; j < n_members; ++j)
        for(size_t k = 0; k < n_members; ++k)
          sigma(k, j) += var_ad[i] * cor_mats[i](k, j);

    Type log_det_sigma;
    matrix<Type> const sigma_inv = atomic::matinvpd(sigma, log_det_sigma);
    Type const entrop_arg = quad_form_sym(rho, lambda);

    term += (
      atomic::logdet(lambda) - quad_form_sym(mu, sigma_inv)
      - mat_mult_trace(lambda, sigma_inv) - log_det_sigma
      + Type(n_members)) / two;
      term -= sqrt_2_pi * quad_form(mu, sigma_inv, delta) + type_M_LN2
        + SNVA::entropy_term(entrop_arg, n_nodes);

    vector<Type> result(1L);
    result[0] = -term;

    // stop recording and return
    ADFun<double> func_out;
    func_out.Dependent(par, result);
    func_out.optimize();
    return func_out;
  })()) { }

  pedigree_element_func_snva(pedigree_element_func_snva&& other):
  n_global(other.n_global), n_members(other.n_members),
  n_priv(other.n_priv),
  this_ad_func(std::move(other.this_ad_func)){ }

  pedigree_element_func_snva(pedigree_element_func_snva const &other):
  n_global(other.n_global), n_members(other.n_members),
  n_priv(other.n_priv) {
    // we may not use the copy constructor but we may use the assignment
    // constructor
    this_ad_func = other.this_ad_func;
  }

  size_t global_dim() const {
    return n_global;
  }
  size_t private_dim() const {
    return n_priv;
  }

  double func(double const *point) const {
    // copy
    // TODO: avoid repeated memory allocation?
    vector<double> par(n_global + n_priv);
    double *p = &par[0];
    for(size_t i = 0; i < n_global + n_priv; ++i, ++p, ++point)
      *p = *point;

    // evaluate
    return this_ad_func.Forward(0, par)[0];
  }
  double grad
  (double const * __restrict__ point, double * __restrict__ gr) const {
    // copy
    // TODO: avoid repeated memory allocation?
    vector<double> par(n_global + n_priv);
    double *p = &par[0];
    for(size_t i = 0; i < n_global + n_priv; ++i, ++p, ++point)
      *p = *point;

    // evaluate
    double const out = this_ad_func.Forward(0, par)[0];
    vector<double> w(1);
    w[0] = 1;

    vector<double> grad = this_ad_func.Reverse(1, w);
    for(size_t i = 0; i < n_global + n_priv; ++i, ++gr)
      *gr = grad[i];

    return out;
  }

  bool thread_safe() const {
    return true;
  }
};

template<class outT> size_t get_n_va(size_t const n_members){
  return 0L;
}

template<> size_t get_n_va<pedigree_element_func_snva>
(size_t const n_members){
  return 2L * n_members + (n_members * (n_members + 1L)) / 2L;
}

template<> size_t get_n_va<pedigree_element_func_gva>
(size_t const n_members){
  return n_members + (n_members * (n_members + 1L)) / 2L;
}
} // namespaces

template <class outT>
SEXP get_pedigree_funcs_generic
  (Rcpp::List data, int const n_nodes, std::string const &link,
   arma::vec const &omega, arma::vec const &beta, arma::vec const &log_sds,
   arma::vec const &va_par, double const eps, double const kappa,
   unsigned const n_threads){
  using ret_T =
    PSQN::optimizer<outT, PSQN::R_reporter, PSQN::R_interrupter>;
  shut_up();
  setup_parallel_ad setup_ADd(n_threads);

  std::vector<outT> funcs;
  size_t const n_clusters = data.size();
  funcs.reserve(n_clusters);

  size_t va_i(0);
  for(size_t g = 0; g < n_clusters; ++g){
    Rcpp::List data_g = data[g];
    size_t const n_members = Rcpp::NumericVector(
      static_cast<SEXP>(data_g["event"])).size(),
                 n_va      = get_n_va<outT>(n_members);
    arma::vec va_par_g(n_va);
    for(size_t i = 0; i < n_va; ++i, ++va_i)
      va_par_g[i] = va_par[va_i];
    funcs.emplace_back(
      data_g, n_nodes, link, omega, beta, log_sds,
      va_par_g, eps, kappa);
  }

  return Rcpp::XPtr<ret_T>(new ret_T(funcs, n_threads));
}

// [[Rcpp::export(rng = false)]]
SEXP get_pedigree_funcs
  (Rcpp::List data, int const n_nodes, std::string const &link,
   arma::vec const &omega, arma::vec const &beta, arma::vec const &log_sds,
   arma::vec const &va_par, double const eps, double const kappa,
   unsigned const n_threads, std::string const &method){
  if(method == "SNVA"){
    return get_pedigree_funcs_generic<pedigree_element_func_snva>
    (data, n_nodes, link, omega, beta, log_sds, va_par, eps, kappa,
     n_threads);
  } else if(method != "GVA")
    throw std::invalid_argument("get_pedigree_funcs: unkown method");

  return get_pedigree_funcs_generic<pedigree_element_func_gva>
    (data, n_nodes, link, omega, beta, log_sds, va_par, eps, kappa,
     n_threads);
}

template <class outT>
Rcpp::List psqn_optim_pedigree_generic
  (Rcpp::NumericVector val, SEXP ptr, double const rel_eps,
   unsigned const max_it, unsigned const n_threads, double const c1,
   double const c2, bool const use_bfgs, int const trace,
   double const cg_tol, bool const strong_wolfe, size_t const max_cg,
   int const pre_method){
  using ret_T =
    PSQN::optimizer<outT, PSQN::R_reporter, PSQN::R_interrupter>;
  Rcpp::XPtr<ret_T> optim(ptr);

  // TODO: check the size

  Rcpp::NumericVector par = clone(val);
  optim->set_n_threads(n_threads);
  auto res = optim->optim(&par[0], rel_eps, max_it, c1, c2,
                          use_bfgs, trace, cg_tol, strong_wolfe, max_cg,
                          static_cast<PSQN::precondition>(pre_method));
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
Rcpp::List psqn_optim_pedigree
  (Rcpp::NumericVector val, SEXP ptr, double const rel_eps,
   unsigned const max_it, unsigned const n_threads, double const c1,
   double const c2, bool const use_bfgs, int const trace,
   double const cg_tol, bool const strong_wolfe, std::string const &method,
   size_t const max_cg, int const pre_method){
  if(method == "SNVA"){
    return psqn_optim_pedigree_generic<pedigree_element_func_snva>
    (val, ptr, rel_eps, max_it, n_threads, c1, c2, use_bfgs, trace, cg_tol,
     strong_wolfe, max_cg, pre_method);
  } else if(method != "GVA")
    throw std::invalid_argument("psqn_optim_pedigree: unkown method");

  return psqn_optim_pedigree_generic<pedigree_element_func_gva>
    (val, ptr, rel_eps, max_it, n_threads, c1, c2, use_bfgs, trace, cg_tol,
     strong_wolfe, max_cg, pre_method);
}

template <class outT>
Rcpp::NumericVector psqn_optim_pedigree_private_generic
  (Rcpp::NumericVector val, SEXP ptr, double const rel_eps,
   unsigned const max_it, unsigned const n_threads, double const c1,
   double const c2){
  using ret_T =
    PSQN::optimizer<outT, PSQN::R_reporter, PSQN::R_interrupter>;
  Rcpp::XPtr<ret_T> optim(ptr);

  // TODO: check the size

  Rcpp::NumericVector par = clone(val);
  optim->set_n_threads(n_threads);
  double const res = optim->optim_priv(&par[0], rel_eps, max_it, c1, c2);
  par.attr("value") = res;
  return par;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector psqn_optim_pedigree_private
  (Rcpp::NumericVector val, SEXP ptr, double const rel_eps,
   unsigned const max_it, unsigned const n_threads, double const c1,
   double const c2, std::string const &method){
  if(method == "SNVA"){
    return psqn_optim_pedigree_private_generic<pedigree_element_func_snva>
    (val, ptr, rel_eps, max_it, n_threads, c1, c2);
  } else if(method != "GVA")
    throw std::invalid_argument("psqn_optim_pedigree_private: unkown method");

  return psqn_optim_pedigree_private_generic<pedigree_element_func_gva>
    (val, ptr, rel_eps, max_it, n_threads, c1, c2);
}

template <class outT>
double eval_psqn_pedigree_generic
  (Rcpp::NumericVector val, SEXP ptr, unsigned const n_threads){
  using ret_T =
    PSQN::optimizer<outT, PSQN::R_reporter, PSQN::R_interrupter>;
  Rcpp::XPtr<ret_T> optim(ptr);

  // TODO: check the size

  optim->set_n_threads(n_threads);
  return optim->eval(&val[0], nullptr, false);
}

// [[Rcpp::export(rng = false)]]
double eval_psqn_pedigree(Rcpp::NumericVector val, SEXP ptr,
                          unsigned const n_threads,
                          std::string const &method){
  if(method == "SNVA"){
    return eval_psqn_pedigree_generic<pedigree_element_func_snva>
    (val, ptr, n_threads);
  } else if(method != "GVA")
    throw std::invalid_argument("eval_psqn_pedigree: unkown method");

  return eval_psqn_pedigree_generic<pedigree_element_func_gva>
    (val, ptr, n_threads);
}

template <class outT>
Rcpp::NumericVector grad_psqn_pedigree_generic
  (Rcpp::NumericVector val, SEXP ptr, unsigned const n_threads){
  using ret_T =
    PSQN::optimizer<outT, PSQN::R_reporter, PSQN::R_interrupter>;
  Rcpp::XPtr<ret_T> optim(ptr);

  // TODO: check the size

  Rcpp::NumericVector grad(val.size());
  optim->set_n_threads(n_threads);
  grad.attr("value") = optim->eval(&val[0], &grad[0], true);

  return grad;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector grad_psqn_pedigree
  (Rcpp::NumericVector val, SEXP ptr, unsigned const n_threads,
   std::string const &method){
  if(method == "SNVA"){
    return grad_psqn_pedigree_generic<pedigree_element_func_snva>
    (val, ptr, n_threads);
  } else if(method != "GVA")
    throw std::invalid_argument("grad_psqn_pedigree: unkown method");

  return grad_psqn_pedigree_generic<pedigree_element_func_gva>
    (val, ptr, n_threads);
}
