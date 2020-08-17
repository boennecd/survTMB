#include <string>
#include "gva-utils.h"
#include "snva-utils.h"
#include "utils.h"

namespace {
/* the atomic_base class has a static list here
 *   https://github.com/kaskr/adcomp/blob/07678703541b81966f18f134d7d7cbf5df9e3345/TMB/inst/include/cppad/local/atomic_base.hpp#L59
 *
 * which is written to here in the destructor
 *   https://github.com/kaskr/adcomp/blob/07678703541b81966f18f134d7d7cbf5df9e3345/TMB/inst/include/cppad/local/atomic_base.hpp#L211
 *
 * Because of the way the `get_cached` functions are written, this mean that
 * each type has to be constructed __before__ calling `get_cached` such that
 * the list gets destructed last. See
 *   https://stackoverflow.com/a/335746/5861244
 *
 * TODO: I only think we need to create one atomic class as the list should
 *       be shared between all derived classes.
 */
template<template<class> class Int>
void fix_atomic_seqfault_helper(){
  Int<      double    > have_to      ("avoid seqfault", 2);
  Int<   AD<double  > > have_to_ad   ("avoid seqfault", 2);
  Int<AD<AD<double> > > have_to_ad_ad("avoid seqfault", 2);
}

template<template<class> class Int>
void get_cached_atomic_objs(size_t const n_nodes){
  Int<   double  >::get_cached(n_nodes);
  Int<AD<double> >::get_cached(n_nodes);
}

} // namespace

// [[Rcpp::export(rng = false)]]
void fix_atomic_seqfault(){
  fix_atomic_seqfault_helper<GaussHermite::GVA::mlogit_integral_atomic>();
  fix_atomic_seqfault_helper<GaussHermite::GVA::probit_integral_atomic>();

  fix_atomic_seqfault_helper<GaussHermite::SNVA::entropy_term_integral>();
  fix_atomic_seqfault_helper<GaussHermite::SNVA::mlogit_integral_atomic>();
  fix_atomic_seqfault_helper<GaussHermite::SNVA::probit_integral_atomic>();
}

// [[Rcpp::export(rng = false)]]
void setup_atomic_cache(size_t const n_nodes, std::string const type,
                        std::string const link = ""){
  if(CppAD::thread_alloc::in_parallel())
    throw std::runtime_error("setup_atomic_cache called in parallel mode");

  if(type == "Laplace") {
    return;
  } else if(type == "GVA"){
    if(link == "PO")
      get_cached_atomic_objs
      <GaussHermite::GVA::mlogit_integral_atomic>(n_nodes);
    else if(link == "probit")
      get_cached_atomic_objs
      <GaussHermite::GVA::probit_integral_atomic>(n_nodes);
    else if(link == "PH" or link.empty()) { }
    else
      throw std::invalid_argument("unkown link (GVA)");

    return;

  } else if(type == "SNVA"){
    GaussHermite::SNVA::entropy_term_integral<   double  >
      have_to   ("avoid seqfaul", n_nodes);
    GaussHermite::SNVA::entropy_term_integral<AD<double> >
      have_ad_to("avoid seqfaul", n_nodes);

    GaussHermite::SNVA::entropy_term_integral<   double  >
      ::get_cached(n_nodes);
    GaussHermite::SNVA::entropy_term_integral<AD<double> >
      ::get_cached(n_nodes);

    if(link == "PO")
      get_cached_atomic_objs
      <GaussHermite::SNVA::mlogit_integral_atomic>(n_nodes);
    else if(link == "probit")
      get_cached_atomic_objs
      <GaussHermite::SNVA::probit_integral_atomic>(n_nodes);
    else if(link == "PH" or link.empty()) { }
    else
      throw std::invalid_argument("unkown link (SNVA)");

    return;

  }

  throw std::invalid_argument("unkown type");
}
