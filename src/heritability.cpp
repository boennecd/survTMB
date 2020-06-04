#define INCLUDE_RCPP
#include "tmb_includes.h"
#include "get-x.h"
#include "snva-utils.h"

namespace {
using namespace GaussHermite::SNVA;

template<class Tout, class T>
vector<Tout> get_args(vector<T> const &omega, vector<T> const &beta,
                      vector<T> const &log_sds, vector<T> const &va_par){
  vector<Tout> out(omega.size() + beta.size() +log_sds.size() + va_par.size());
  Tout *o = &out[0];
  for(int i = 0; i < omega.size(); ++i, ++o)
    *o = Tout(omega[i]);
  for(int i = 0; i < beta.size(); ++i, ++o)
    *o = Tout(beta[i]);
  for(int i = 0; i < log_sds.size(); ++i, ++o)
    *o = Tout(log_sds[i]);
  for(int i = 0; i < va_par.size(); ++i, ++o)
    *o = Tout(va_par[i]);

  return out;
}

/* object to hold the data for a given cluster */
template<class Type>
struct cluster_data {
  Rcpp::List data;

  const DATA_VECTOR(y);
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

  size_t const n_members = y.size();

  cluster_data(Rcpp::List data): data(data) {
    if((size_t)event.size() != n_members)
      throw std::invalid_argument("cluster_data<Type>: invalid event");
    else if((size_t)X.cols()  != n_members)
      throw std::invalid_argument("cluster_data<Type>: invalid X");
    else if((size_t)XD.cols() != n_members)
      throw std::invalid_argument("cluster_data<Type>: invalid XD");
    else if((size_t)Z.cols() != n_members)
      throw std::invalid_argument("cluster_data<Type>: invalid Z");
    for(auto &V : cor_mats)
      if((size_t)V.rows() != n_members or (size_t)V.cols() != n_members)
        throw std::invalid_argument("cluster_data<Type>: invalid cor_mats");

    data = Rcpp::List();
  }
};

template<class Type>
class VA_worker {
  Rcpp::List data, parameters;

  std::vector<cluster_data<Type> > const c_data = ([&](){
    Rcpp::List c_data_R = data["c_data"];

    std::vector<cluster_data<Type> > out;
    out.reserve(c_data_R.size());
    for(auto x : c_data_R)
      out.emplace_back(Rcpp::List(x));

    return out;
  })();

  const DATA_INTEGER(n_threads);
  const DATA_LOGICAL(sparse_hess);
  const DATA_INTEGER(n_nodes);
  const DATA_STRING(link);

  const PARAMETER_VECTOR(omega);
  const PARAMETER_VECTOR(beta);
  const PARAMETER_VECTOR(log_sds);
  const PARAMETER_VECTOR(va_par);
  const PARAMETER(eps);
  const PARAMETER(kappa);

public:
  size_t const d_o = omega.size(),
            n_mats = log_sds.size(),
        n_clusters = c_data.size(),
             n_obs = ([&]{
               size_t out(0L);
               for(auto &x : c_data)
                 out += x.y.size();

               return out;
             })(),
            n_pars = omega.size() + beta.size() + log_sds.size() +
              va_par.size();

#ifdef _OPENMP
  std::size_t const n_blocks = n_threads;
#else
  std::size_t const n_blocks = 1L;
#endif

  VA_worker(Rcpp::List data, Rcpp::List parameters):
  data(data), parameters(parameters) {
    /* checks */
    {
      size_t expected_va_pars(0L);
      for(auto &x : c_data)
        expected_va_pars +=
          2L * x.n_members + (x.n_members * (x.n_members + 1L)) / 2L;

      if(n_nodes < 1L)
        throw std::invalid_argument("VA_worker<Type>: invalid n_nodes");
      else if(n_threads < 1L)
        throw std::invalid_argument("VA_worker<Type>: invalid n_threads");
      else if((size_t)va_par.size() != expected_va_pars)
        throw std::invalid_argument("VA_worker<Type>: invalid va_par");
      for(auto &x : c_data)
        if((size_t)x.X.rows() != d_o)
          throw std::invalid_argument("VA_worker<Type>: invalid c_data (X)");
        else if(x.Z.rows() != beta.size())
          throw std::invalid_argument("VA_worker<Type>: invalid c_data (Z)");
        else if(x.cor_mats.size() != n_mats)
          throw std::invalid_argument(
              "VA_worker<Type>: invalid c_data (cor_mats)");
        else if((size_t)x.XD.rows() != d_o)
          throw std::invalid_argument("VA_worker<Type>: invalid c_data (XD)");
    }

#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    data = Rcpp::List();
    parameters = Rcpp::List();
  }

  template<typename Tout>
  vector<Tout> get_args() const {
    return ::get_args<Tout, Type>(omega, beta, log_sds, va_par);
  }

  Type operator()(vector<Type> &args) const {
    /* assign constant */
    Type const sqrt_2_pi(sqrt(M_2_PI)),
                     one(1.),
                     two(2.),
                   small(std::numeric_limits<double>::epsilon()),
              type_M_LN2(M_LN2);

    /* assign parameters */
    if((size_t)args.size() != n_pars)
      throw std::invalid_argument("VA_worker::operator(): invalid args");

    Type const *a = &args[0];
    auto set_vec = [&](size_t const n){
      vector<Type> out(n);
      for(size_t i = 0; i < n; ++i, ++a)
        out[i] = *a;
      return out;
    };

    vector<Type> const aomega = set_vec(omega.size()),
                        abeta = set_vec(beta.size()),
                     alog_sds = set_vec(log_sds.size()),
                        a_var = ([&](){
                         vector<Type> out(alog_sds.size());
                         for(size_t i = 0; i < (size_t)out.size(); ++i)
                           out[i] = exp(two * alog_sds[i]);
                         return out;
                         })();

    /* handle VA pars */
    std::vector<SNVA_MD_input<Type> > ava_par;
    ava_par.reserve(c_data.size());
    for(auto &c_dat : c_data){
      size_t const n_ele = c_dat.n_members,
                n_parmas = 2L * n_ele + (n_ele * (n_ele + 1L)) / 2L;
      ava_par.emplace_back(
        SNVA_MD_theta_DP_to_DP(a, n_parmas, n_ele));
      a += n_parmas;
    }

    survTMB::accumulator_mock<Type> result;
    bool const is_in_parallel = CppAD::thread_alloc::in_parallel();
    ph    <Type>     ph_func(eps, kappa, n_nodes);
    po    <Type>     po_func(eps, kappa, n_nodes);
    probit<Type> probit_func(eps, kappa, n_nodes);
    for(size_t g = 0; g < n_clusters; ++g){
      if(is_in_parallel and !is_my_region(*result.obj)){
        result.obj->parallel_region();
        continue;
      }

      /* compute objects needed for the variational distribution */
      matrix<Type> const &lambda = ava_par[g].va_lambdas[0];
      vector<Type> const &mu = ava_par[g].va_mus[0],
                        &rho = ava_par[g].va_rhos[0];

      vector<Type> const delta = ([&](){
        vector<Type> out = lambda * rho;
        Type const denom = sqrt(one + vec_dot(out, rho));
        out /= denom;
        return out;
      })();

      /* add terms from each cluster member */
      Type term(0.);
      auto const &c_dat = c_data[g];
      size_t const n_members =  c_dat.n_members;
      for(size_t i = 0; i < n_members; ++i){
        vector<Type> const x  = c_dat.X .col(i),
                           xd = c_dat.XD.col(i),
                           z  = c_dat.Z .col(i);

        Type const &sd_sq = lambda(i, i),
                       sd = sqrt(sd_sq),
                       &d = delta(i),
                      rho = d / sd_sq / sqrt(one - d * d / sd_sq),
                 d_scaled = sqrt_2_pi * d,
                dist_mean = mu(i) + d_scaled,
                 dist_var = sd_sq - d_scaled * d_scaled;

        if(link == "PH")
          term += ph_func(
            vec_dot(aomega, x) + vec_dot(abeta, z), vec_dot(aomega, xd),
            c_dat.event[i], mu[i], sd, rho, d, sd_sq, dist_mean, dist_var);
        else if(link == "PO")
          term += po_func(
            vec_dot(aomega, x) + vec_dot(abeta, z), vec_dot(aomega, xd),
            c_dat.event[i], mu[i], sd, rho, d, sd_sq, dist_mean, dist_var);
        else if(link == "probit")
          term += probit_func(
            vec_dot(aomega, x) + vec_dot(abeta, z), vec_dot(aomega, xd),
            c_dat.event[i], mu[i], sd, rho, d, sd_sq, dist_mean, dist_var);
        else
          error("'%s' not implemented", link.c_str());
      }

      /* add prior and entropy terms */
      matrix<Type> sigma(n_members, n_members);
      sigma.setZero();
      for(int i = 0; i < a_var.size(); ++i)
        sigma += a_var[i] * c_dat.cor_mats[i];

      Type log_det_sigma;
      matrix<Type> const sigma_inv = atomic::matinvpd(sigma, log_det_sigma);

      term += (
        atomic::logdet(lambda) - quad_form_sym(mu, sigma_inv) -
          mat_mult_trace(lambda, sigma_inv) - log_det_sigma +
          Type(n_members)) / two;
      term -= sqrt_2_pi * quad_form(mu, sigma_inv, delta) + type_M_LN2;

      result -= term;
    }

    return result;
  }
};

class VA_func {
  using ADd   = CppAD::AD<double>;
  using ADdd  = CppAD::AD<ADd>;
  using ADddd = CppAD::AD<ADdd>;
  template<class Type>
  using ADFun = CppAD::ADFun<Type>;

  size_t n_pars;

public:
  size_t get_n_pars() const {
    return n_pars;
  }

  std::vector<std::unique_ptr<ADFun<double> > >   funcs;

  VA_func(Rcpp::List data, Rcpp::List parameters){
    {
      /* to compute function and gradient */
      VA_worker<ADd> w(data, parameters);
      funcs.resize(w.n_blocks);
      n_pars = w.n_pars;
      vector<ADd> args = w.get_args<ADd>();

#ifdef _OPENMP
#pragma omp parallel for if(w.n_blocks > 1L) firstprivate(args)
#endif
      for(unsigned i = 0; i < w.n_blocks; ++i){
        funcs[i].reset(new ADFun<double>());

        CppAD::Independent(args);
        vector<ADd> y(1);
        y[0] = w(args);

        funcs[i]->Dependent(args, y);
        funcs[i]->optimize();
      }
    }
  }
};

} // namespaces

// [[Rcpp::export(rng = false)]]
SEXP get_herita_funcs
  (Rcpp::List data, Rcpp::List parameters){
  shut_up();

  unsigned const n_threads(data["n_threads"]);
  setup_parallel_ad setup_ADd(n_threads);
  return Rcpp::XPtr<VA_func>(new VA_func(data, parameters));
}

// [[Rcpp::export(rng = false)]]
double herita_funcs_eval_lb(SEXP p, SEXP par){
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  std::vector<std::unique_ptr<CppAD::ADFun<double> > > &funcs = ptr->funcs;
  vector<double> parv = get_vec<double>(par);
  if((size_t)parv.size() != ptr->get_n_pars())
    throw std::invalid_argument("herita_funcs_eval_lb: invalid par");

  unsigned const n_blocks = ptr->funcs.size();
  double out(0);
#ifdef _OPENMP
  /* #pragma omp parallel for if(n_blocks > 1L) firstprivate(parv) reduction(+:out) */
#pragma omp parallel for ordered if(n_blocks > 1L) firstprivate(parv)  schedule(static, 1)
#endif
  for(unsigned i = 0; i < n_blocks; ++i){
    double const term = funcs[i]->Forward(0, parv)[0];
#ifdef _OPENMP
    /* TODO: replace with a reduction */
#pragma omp ordered
#endif
    out += term;
  }

  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector herita_funcs_eval_grad(SEXP p, SEXP par){
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  std::vector<std::unique_ptr<CppAD::ADFun<double> > > &funcs = ptr->funcs;
  vector<double> parv = get_vec<double>(par);
  if((size_t)parv.size() != ptr->get_n_pars())
    throw std::invalid_argument("herita_funcs_eval_grad: invalid par");

  unsigned const n_blocks = ptr->funcs.size();
  vector<double> grad(parv.size());
  grad.setZero();

#ifdef _OPENMP
#pragma omp parallel for ordered if(n_blocks > 1L) firstprivate(parv) schedule(static, 1)
#endif
  for(unsigned i = 0; i < n_blocks; ++i){
    funcs[i]->Forward(0, parv);
    vector<double> w(1);
    w[0] = 1;

    vector<double> grad_i = funcs[i]->Reverse(1, w);
#ifdef _OPENMP
    /* TODO: replace with a reduction */
#pragma omp ordered
#endif
    grad += grad_i;
  }

  std::size_t const n = grad.size();
  Rcpp::NumericVector out(n);
  for(unsigned i = 0; i < n; ++i)
    out[i] = grad[i];

  return out;
}
