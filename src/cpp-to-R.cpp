#define INCLUDE_RCPP
#include "tmb_includes.h"
#include "fastgl.h"
#include "get-x.h"
#include "clear-mem.h"
#ifdef _OPENMP
#include "omp.h"
#endif

// [[Rcpp::export(rng = false)]]
Rcpp::List get_gl_rule(unsigned const n){
  if(n == 0L)
    throw std::invalid_argument("get_gl_rule: n is zero");

  auto const &dat = fastgl::GLPairsCached<double>(n);
  Rcpp::NumericVector x(n), w(n);
  for(unsigned i = 0; i < n; ++i){
    auto const &dat_i = dat[i];
    x[i] = dat_i.x;
    w[i] = dat_i.weight;
  }

  return Rcpp::List::create(
    Rcpp::Named("node")   = x,
    Rcpp::Named("weight") = w);
}


//' Clears the Memory used by CppAD
//' Clear all pointers to CppAD objects and deallocate all memory used by
//' CppAD. This will make all subsequent use C++ pointer invalid.
//'
//' @param max_n_threads maximum number of threads which have been used.
//' @param keep_work_space logical for whether to keep
//'
//' @return An integer which is one if all memory could be freed.
//'
//' @export
// [[Rcpp::export(rng = false)]]
int clear_cppad_mem
  (unsigned const max_n_threads = 1L, bool const keep_work_space = false){
  clear_clearables();

  if(keep_work_space)
    return CppAD::thread_alloc::free_all();

  {
    setup_parallel_ad setup_ADd(max_n_threads);
    for(size_t i = 0; i < max_n_threads; ++i)
      CppAD::thread_alloc::free_available(i);

    add_atomics_to_be_cleared();
    clear_atomics<         double>      ();
    clear_atomics<      AD<double> >    ();
    clear_atomics<   AD<AD<double> > >  ();
    clear_atomics<AD<AD<AD<double> > > >();
  }

  setup_parallel_ad setup_ADd(1L);
  return CppAD::thread_alloc::free_all();
}

// [[Rcpp::export(rng = false)]]
int set_n_threads(int const n_threads){
#ifdef _OPENMP
  if(n_threads > 0L)
    omp_set_num_threads(n_threads);
  return omp_get_max_threads();
#else
  return 1L;
#endif
}
