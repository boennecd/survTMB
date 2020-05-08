#include <Rcpp.h>
#include "fastgl.h"

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
