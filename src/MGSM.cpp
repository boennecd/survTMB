#define INCLUDE_RCPP
#include "get-x.h"
#include "utils.h"
#include "snva.h"
#include "gva.h"
#include <memory>
#include <vector>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

template<typename Tout>
vector<Tout> get_args_va(Rcpp::List parameters) {
  using Type = double;
  PARAMETER(eps);
  PARAMETER(kappa);
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(theta_va);

  return get_args_va<Tout, Type>(eps, kappa, b, theta, theta_va);
}

/* define and declare functor to use */
template<class Type>
class VA_worker {
  Rcpp::List data, parameters;

  SETUP_DATA;
  PARAMETER_VECTOR(theta_VA);
  DATA_INTEGER(n_nodes);
  DATA_STRING(param_type);

  std::size_t const n_b = b       .size(),
                    n_t = theta   .size(),
                    n_v = theta_VA.size(),
                 n_para = 2L + n_b + n_t + n_v;

public:
#ifdef _OPENMP
  std::size_t const n_blocks = n_threads;
#else
  std::size_t const n_blocks = 1L;
#endif

  VA_worker(Rcpp::List data, Rcpp::List parameters):
  data(data), parameters(parameters) {
    SETUP_DATA_CHECK;

#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif

    data = Rcpp::List();
    parameters = Rcpp::List();
  }

  template<typename Tout>
  vector<Tout> get_args_va() const {
    return ::get_args_va<Tout, Type>(eps, kappa, b, theta, theta_VA);
  }

  Type operator()(vector<Type> &args) const {
    if((unsigned)args.size() != n_para)
      error("VA_worker: invalid args length");
    Type eps = args[0],
       kappa = args[1];
    Type const *a = &args[2];

    vector<Type> b(n_b);
    for(unsigned i = 0; i < n_b; ++i, ++a)
      b[i] = *a;
    vector<Type> theta(n_t);
    for(unsigned i = 0; i < n_t; ++i, ++a)
      theta[i] = *a;
    vector<Type> theta_VA(n_v);
    for(unsigned i = 0; i < n_v; ++i, ++a)
      theta_VA[i] = *a;

    survTMB::accumulator_mock<Type> result;
    if(app_type =="GVA"){
      GVA(COMMON_CALL, theta_VA, n_nodes);
      return result;

    } else if(app_type == "SNVA"){
      SNVA(COMMON_CALL, theta_VA, n_nodes, param_type);
      return result;

    }

    error("VA_worker<Type>: approximation method '%s' is not implemented",
          app_type.c_str());
    return Type(0);
  }
};

struct VA_func {
  using ADd   = CppAD::AD<double>;
  using ADdd  = CppAD::AD<ADd>;
  using ADddd = CppAD::AD<ADdd>;
  template<class Type>
  using ADFun = CppAD::ADFun<Type>;

                  std::vector<std::unique_ptr<ADFun<double> > >   funcs;
  std::unique_ptr<std::vector<std::unique_ptr<ADFun<double> > > > grads;

  struct sparse_mat_data {
    vector<int> row_idx, col_idx;
    CppAD::ADFun<double> ddf;
  };
  std::unique_ptr<sparse_mat_data> sparse_hess_dat;

  VA_func(Rcpp::List data, Rcpp::List parameters){
    {
      /* to compute function and gradient */
      VA_worker<ADd> w(data, parameters);
      funcs.resize(w.n_blocks);
      vector<ADd> args = w.get_args_va<ADd>();

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

    /* TODO: build this on request afterwards */
    DATA_LOGICAL(dense_hess);
    if(dense_hess){
      /* to compute dense Hessian
       * TODO: use base2ad if CppAD gets updated */
      grads.reset(new std::vector<std::unique_ptr<ADFun<double> > >());
      auto &grs = *grads;

      VA_worker<ADdd> w(data, parameters);
      grs.resize(w.n_blocks);
      vector<ADdd> x = w.get_args_va<ADdd>();
#ifdef _OPENMP
#pragma omp parallel for if(w.n_blocks > 1L) firstprivate(x)
#endif
      for(unsigned i = 0; i < w.n_blocks; ++i){
        grs[i].reset(new ADFun<double>());

        CppAD::Independent(x);
        vector<ADdd> y(1);
        y[0] = w(x);

        ADFun<ADd> tmp;
        tmp.Dependent(x, y);
        tmp.optimize();

        vector<ADd> xx(x.size());
        for(unsigned i = 0; i < x.size(); ++i)
          xx[i] = CppAD::Value(x[i]);

        CppAD::Independent(xx);
        vector<ADd> yy = tmp.Jacobian(xx);

        grs[i]->Dependent(xx, yy);
        grs[i]->optimize();
      }
    }

    /* TODO: build this on request afterwards */
    DATA_LOGICAL(sparse_hess);
    if(sparse_hess){
      /* to compute sparse Hessian
       * TODO: use subgraph_jac_rev if CppAD gets updated */
      VA_worker<ADddd> w(data, parameters);
      vector<ADddd> x = w.get_args_va<ADddd>();
      unsigned const n_vars = x.size();

      /* record f */
      CppAD::Independent(x);
      vector<ADddd> y(1);
      y[0] = w(x);
      CppAD::ADFun<ADdd> f;
      f.Dependent(x, y);
      f.optimize();

      /* record f'' */
      vector<ADdd> xx(n_vars),
                   wi(1);
      for(unsigned i = 0; i < n_vars; ++i)
        xx[i] = CppAD::Value(x[i]);
      CppAD::Independent(xx);
      f.Forward(0, xx);
      wi[0] = 1.;
      vector<ADdd> yy = f.Reverse(1, wi);

      CppAD::ADFun<ADd> df;
      df.Dependent(xx, yy);
      df.optimize();

      /* record f'' but in a sparse manner. Thus, we first need to
       * find the number of non-zero entries and their positions */
      vector<bool> keepcol(n_vars);
      for(unsigned i = 0; i < 2; ++i)
        keepcol[i] = false;
      for(unsigned i = 2; i < n_vars; ++i)
        keepcol[i] = true;
      df.my_init(keepcol);

      auto keep_col = [&](unsigned const col){
        return keepcol[col];
      };
      auto keep_row = [&](unsigned const row, unsigned const col){
        return keep_col(col) and row >= col;
      };

      unsigned colisize,
               n_non_zero(0); // Count number of non-zeros (m)
      for(unsigned i = 0; i < df.colpattern.size(); i++){
        colisize = df.colpattern[i].size();
        if(keep_col(i))
          for(unsigned j = 0; j < colisize; j++)
            n_non_zero += keep_row(df.colpattern[i][j], i);
      }

      // Allocate index vectors of non-zero pairs
      vector<int> row_idx(n_non_zero);
      vector<int> col_idx(n_non_zero);
      // Prepare reverse sweep for Hessian columns
      vector<ADd> u(n_vars);
      vector<ADd> v(n_vars);
      for(unsigned i = 0; i < n_vars; i++)
        v[i] = 0.0;
      vector<ADd> xxx(n_vars);
      for(unsigned i = 0; i < n_vars; i++)
        xxx[i]  = CppAD::Value(CppAD::Value(xx[i]));
      vector<ADd> yyy(n_non_zero);

      // Do sweeps and fill in non-zero index pairs
      CppAD::Independent(xxx);
      df.Forward(0, xxx);

      unsigned idx(0);
      for(unsigned i = 0; i < n_vars; i++){
        if(keep_col(i)) {
          df.myReverse(1, v, i /*range comp*/, u /*domain*/);
          CppAD::vector<int> &icol = df.colpattern[i];

          for(unsigned j = 0; j < icol.size(); j++){
            if(keep_row(icol[j], i)){
              row_idx[idx] = icol[j];
              col_idx[idx] = i;
              yyy    [idx] = u[icol[j]];
              idx++;
            }
          }
        }
      }

      /* store output */
      sparse_hess_dat.reset(new sparse_mat_data());
      auto &shd = *sparse_hess_dat;
      shd.ddf.Dependent(xxx, yyy);
      shd.row_idx = std::move(row_idx);
      shd.col_idx = std::move(col_idx);
    }
  }
};

} // namespace

// [[Rcpp::export(rng = false)]]
SEXP get_VA_funcs
  (Rcpp::List data, Rcpp::List parameters){
  shut_up();

  unsigned const n_threads(data["n_threads"]);
  setup_parallel_ad<double>
    setup_ADd (n_threads);
  setup_parallel_ad<CppAD::AD<double> >
    setup_ADdd(n_threads, false);

  return Rcpp::XPtr<VA_func>(new VA_func(data, parameters));
}

// [[Rcpp::export(rng = false)]]
double VA_funcs_eval_lb
  (SEXP p, SEXP par){
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  std::vector<std::unique_ptr<CppAD::ADFun<double> > > &funcs = ptr->funcs;
  vector<double> parv = get_vec<double>(par);

  unsigned const n_blocks = ptr->funcs.size();
  double out(0);
#ifdef _OPENMP
#pragma omp parallel for if(n_blocks > 1L) firstprivate(parv) reduction(+:out)
#endif
  for(unsigned i = 0; i < n_blocks; ++i)
    out += funcs[i]->Forward(0, parv)[0];

  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector VA_funcs_eval_grad
  (SEXP p, SEXP par){
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  std::vector<std::unique_ptr<CppAD::ADFun<double> > > &funcs = ptr->funcs;
  vector<double> parv = get_vec<double>(par);

  unsigned const n_blocks = ptr->funcs.size();
  vector<double> grad(parv.size());
  grad.setZero();

#ifdef _OPENMP
#pragma omp parallel for if(n_blocks > 1L) firstprivate(parv)
#endif
  for(unsigned i = 0; i < n_blocks; ++i){
    funcs[i]->Forward(0, parv);
    vector<double> w(1);
    w[0] = 1;

    vector<double> grad_i = funcs[i]->Reverse(1, w);
#ifdef _OPENMP
    /* TODO: replace with a reduction */
#pragma omp critical
#endif
    grad += grad_i;
  }

  std::size_t const n = grad.size();
  Rcpp::NumericVector out(n);
  for(unsigned i = 0; i < n; ++i)
    out[i] = grad[i];

  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix VA_funcs_eval_hess
  (SEXP p, SEXP par){
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  if(!ptr->grads)
    throw std::invalid_argument(
        "VA_funcs_eval_hess: no ADFun objects to compute the Hessian");

  std::vector<std::unique_ptr<CppAD::ADFun<double> > >
    &grads = *ptr->grads;

  vector<double> parv = get_vec<double>(par);
  unsigned const n_blocks = grads.size(),
                 n_vars   = parv.size();
  vector<double> hess(n_vars * n_vars);
  hess.setZero();

#ifdef _OPENMP
#pragma omp parallel for if(n_blocks > 1L) firstprivate(parv)
#endif
  for(unsigned i = 0; i < n_blocks; ++i){
    vector<double> hess_i = grads[i]->Jacobian(parv);
#ifdef _OPENMP
    /* TODO: replace with a reduction */
#pragma omp critical
#endif
    hess += hess_i;
  }

  Rcpp::NumericMatrix out(n_vars, n_vars);
  double       *o = &out[0];
  double const *g = &hess[0];
  for(unsigned j = 0; j < n_vars; ++j)
    for(unsigned i = 0; i < n_vars; ++i, ++o, ++g)
      *o = *g;

  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List VA_funcs_eval_hess_sparse
  (SEXP p, SEXP par){
  using Rcpp::Named;
  shut_up();

  Rcpp::XPtr<VA_func > ptr(p);
  if(!ptr->sparse_hess_dat)
    throw std::invalid_argument(
        "VA_funcs_eval_hess_sparse: no object to compute the sparse Hessian");

  vector<double> parv = get_vec<double>(par);
  auto &shd = *ptr->sparse_hess_dat;
  vector<double> val = shd.ddf.Forward(0, parv);

  auto get_integer_vec = [&](vector<int> x){
    unsigned const n = x.size();
    Rcpp::IntegerVector out(n);
    for(unsigned i = 0; i < n; ++i)
      out[i] = x[i];
    return out;
  };
  auto get_numeric_vec = [&](vector<double> x){
    unsigned const n = x.size();
    Rcpp::NumericVector out(n);
    for(unsigned i = 0; i < n; ++i)
      out[i] = x[i];
    return out;
  };

  return Rcpp::List::create(
    Named("val")     = get_numeric_vec(val),
    Named("row_idx") = get_integer_vec(shd.row_idx),
    Named("col_idx") = get_integer_vec(shd.col_idx));
}
