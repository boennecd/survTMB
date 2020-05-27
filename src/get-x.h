#ifndef GET_X_H
#define GET_X_H
#include "tmb_includes.h"

#ifdef _OPENMP
#include <omp.h>
#endif

inline void shut_up(){
  config.trace.atomic = false;
  config.trace.parallel = false;
  config.trace.optimize = false;
}

/* redefine the macros used by TMB */
template<class Type>
vector<Type> get_vec(SEXP obj){
  Rcpp::NumericVector org(obj);
  size_t const n = org.size();
  vector<Type> out(n);

  double const *o = &org(0);
  for(unsigned i = 0; i < n; i++, o++)
    out[i] = Type(*o);

  return out;
}

template<class Type>
matrix<Type> get_mat(SEXP obj){
  Rcpp::NumericMatrix org(obj);
  std::size_t const n = org.nrow(),
                    m = org.ncol();
  matrix<Type> out(n, m);

  double const *o = &org(0, 0);
  for(unsigned j = 0; j < m; j++)
    for(unsigned i = 0; i < n; i++, o++)
      out(i, j) = Type(*o);

  return out;
}

#ifdef DATA_STRING
#undef DATA_STRING
#endif
#define DATA_STRING(name)                                      \
  std::string name = data["" #name ""];                        \

#define DATA_LOGICAL(name)                                     \
  bool name = data["" #name ""];

#ifdef DATA_VECTOR
#undef DATA_VECTOR
#endif
#define DATA_VECTOR(name)                                      \
  vector<Type> name = get_vec<Type>(data["" #name ""]);

#ifdef DATA_MATRIX
#undef DATA_MATRIX
#endif
#define DATA_MATRIX(name)                                      \
  matrix<Type> name = get_mat<Type>(data["" #name ""]);

#ifdef DATA_IVECTOR
#undef DATA_IVECTOR
#endif
#define DATA_IVECTOR(name)                                     \
  vector<int> name = ([&](){                                   \
    Rcpp::IntegerVector org = data["" #name ""];               \
    std::size_t const n = org.size();                          \
    vector<int> out(n);                                        \
                                                               \
    int const *o = &org(0);                                    \
    for(unsigned i = 0; i < n; ++o, ++i)                       \
      out[i] = org(i);                                         \
                                                               \
    return out;                                                \
  })();

#ifdef DATA_INTEGER
#undef DATA_INTEGER
#endif
#define DATA_INTEGER(name)                                     \
  int name = data["" #name ""];

#ifdef PARAMETER
#undef PARAMETER
#endif
#define PARAMETER(name)                                        \
  double name = parameters["" #name ""];

#ifdef PARAMETER_VECTOR
#undef PARAMETER_VECTOR
#endif
#define PARAMETER_VECTOR(name)                                 \
  vector<Type> name =                                          \
    get_vec<Type>(parameters["" #name ""]);                     \

#ifdef PARAMETER_MATRIX
#undef PARAMETER_MATRIX
#endif
#define PARAMETER_MATRIX(name)                                 \
  matrix<Type> name =                                          \
    get_mat<Type>(parameters["" #name ""]);

template<typename Tout, typename Tin>
vector<Tout> get_args_va
  (Tin const &eps, Tin const &kappa, vector<Tin> const &b,
   vector<Tin> const &theta, vector<Tin> const &theta_VA) {
  std::size_t const n_b = b       .size(),
                    n_t = theta   .size(),
                    n_v = theta_VA.size();
  vector<Tout> out(2L + n_b + n_t + n_v);
  out[0] = eps;
  out[1] = kappa;
  Tout *o = &out[2];

  for(unsigned i = 0; i < n_b; ++i, ++o)
    *o = Tout(b[i]);
  for(unsigned i = 0; i < n_t; ++i, ++o)
    *o = Tout(theta[i]);
  for(unsigned i = 0; i < n_v; ++i, ++o)
    *o = Tout(theta_VA[i]);

  return out;
}

#ifdef _OPENMP
inline bool is_in_parallel(){
  return static_cast<bool>(omp_in_parallel());
}
inline size_t get_thread_num(){
  return static_cast<size_t>(omp_get_thread_num());
}
#endif

template<class Type>
struct setup_parallel_ad {
#ifdef _OPENMP
  setup_parallel_ad(std::size_t const nthreads, bool const setup = true) {
    if(setup){
      CppAD::thread_alloc::parallel_setup(
        nthreads, is_in_parallel, get_thread_num);
      CppAD::thread_alloc::hold_memory(true);
    }

    CppAD::parallel_ad<Type>();
  }
#else
  setup_parallel_ad(unsigned const nthreads, bool const setup = true) { }
#endif
};

#endif
