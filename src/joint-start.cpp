#define INCLUDE_RCPP
#include "tmb_includes.h"
#include "splines.h"
#include "fastgl.h"

template<bool grad>
arma::vec joint_start_baseline_inner
  (arma::vec const &Y, arma::vec const &tstart, arma::vec const &tstop,
   arma::vec const &omega, arma::vec const &offsets, unsigned const n_nodes,
   arma::vec const &bound_knots, arma::vec const &inter_knots){
  splines::ns basis(bound_knots, inter_knots, true);
  size_t const n_out = grad ? basis.get_n_basis() : 1L;
  arma::vec out(n_out, arma::fill::zeros);
  auto const &xw = fastgl::GLPairsCached<double>(n_nodes);
  arma::vec wrk(basis.get_n_basis()), term_v(n_out);

  size_t const n = Y.n_elem;
#ifdef DO_CHECKS
  if(tstart.n_elem != n)
    throw std::invalid_argument("joint_start_baseline: invalid tstart");
  if(tstop.n_elem != n)
    throw std::invalid_argument("joint_start_baseline: invalid tstop");
  if(offsets.n_elem != n)
    throw std::invalid_argument("joint_start_baseline: invalid offsets");
  if(bound_knots.n_elem != 2L)
    throw std::invalid_argument("joint_start_baseline: invalid bound_knots");
  if(omega.n_elem != basis.get_n_basis())
    throw std::invalid_argument("joint_start_baseline: invalid omega");
  if(n_nodes == 0L)
    throw std::invalid_argument("joint_start_baseline: invalid n_nodes");
#endif

  for(size_t i = 0; i < n; ++i){
    double term_s(0.);
    if(grad)
      term_v.zeros();
    double const d1 = (tstop[i] - tstart[i]) / 2.,
                 d2 = (tstop[i] + tstart[i]) / 2.;
    for(auto const &xwi : xw){
      double const node = d1 * xwi.x + d2;
      basis(wrk, log(node));
      double const g = xwi.weight * exp(arma::dot(wrk, omega));
      if(grad)
        term_v -= g * wrk;
      else
        term_s -= g;
    }

    if(grad)
      term_v *= d1 * exp(offsets[i]);
    else
      term_s *= d1 * exp(offsets[i]);

    if(Y[i] > 0){
      basis(wrk, log(tstop[i]));
      if(grad)
        term_v += wrk;
      else
        term_s += offsets[i] + arma::dot(wrk, omega);
    }

    if(grad)
      out += term_v;
    else
      out += term_s;
  }

  return out;
}

/**
  Args:
    Y: zero/one vector with event indicators.
    tstart: left truncation time.
    tstop: right-censoring time or event time.
    omega: coefficients.
    offsets: offsets.
    n_nodes: integer with number of Gauss-Legendre quadrature nodes.
    bound_knots: boundary knots.
    inter_knots: interior knots.
    grad: logical for whether to compute the gradient og the log-likelihood.
 */

// [[Rcpp::export(rng = false)]]
arma::vec joint_start_baseline
  (arma::vec const &Y, arma::vec const &tstart, arma::vec const &tstop,
   arma::vec const &omega, arma::vec const &offsets, unsigned const n_nodes,
   arma::vec const &bound_knots, arma::vec const &inter_knots,
   bool const grad){
  if(grad)
    return(joint_start_baseline_inner<true>
             (Y, tstart, tstop, omega, offsets, n_nodes, bound_knots,
              inter_knots));

  return(joint_start_baseline_inner<false>
           (Y, tstart, tstop, omega, offsets, n_nodes, bound_knots,
            inter_knots));
}
