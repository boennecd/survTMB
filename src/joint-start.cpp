#define INCLUDE_RCPP
#include "tmb_includes.h"
#include "splines.h"
#include "fastgl.h"

template<bool grad>
arma::vec joint_start_ll_inner
  (arma::vec const &Y, arma::vec const &tstart, arma::vec const &tstop,
   arma::vec const &omega, arma::vec const &delta, arma::mat const &Z,
   unsigned const n_nodes, arma::vec const &bound_knots,
   arma::vec const &inter_knots){
  splines::ns basis(bound_knots, inter_knots, true);
  size_t const n_out = grad ? basis.get_n_basis() + delta.size() : 1L;
  arma::vec out(n_out, arma::fill::zeros);
  auto const &xw = fastgl::GLPairsCached<double>(n_nodes);
  arma::vec wrk(basis.get_n_basis()), term_v(n_out);

  bool const has_delta = delta.n_elem > 0;
  size_t const o_idx_start = 0L,
                 o_idx_end = omega.n_elem - 1L,
               g_idx_start = o_idx_end + 1L,
                 g_idx_end = o_idx_end + delta.n_elem,
                         n = Y.n_elem;

#ifdef DO_CHECKS
  if(tstart.n_elem != n)
    throw std::invalid_argument("joint_start_ll: invalid tstart");
  else if(tstop.n_elem != n)
    throw std::invalid_argument("joint_start_ll: invalid tstop");
  else if(Z.n_rows != delta.n_elem)
    throw std::invalid_argument("joint_start_ll: invalid Z");
  else if(bound_knots.n_elem != 2L)
    throw std::invalid_argument("joint_start_ll: invalid bound_knots");
  else if(omega.n_elem != basis.get_n_basis() or omega.n_elem < 1L)
    throw std::invalid_argument("joint_start_ll: invalid omega");
  else if(n_nodes == 0L)
    throw std::invalid_argument("joint_start_ll: invalid n_nodes");
#endif

  for(size_t i = 0; i < n; ++i){
    double term_s(0.);
    if(grad)
      term_v.zeros();
    double const d1 = (tstop[i] - tstart[i]) / 2.,
                 d2 = (tstop[i] + tstart[i]) / 2.,
       fixed_effect = has_delta ? arma::dot(Z.col(i), delta) : 0.;
    for(auto const &xwi : xw){
      double const node = d1 * xwi.x + d2;
      basis(wrk, log(node));
      double const g = xwi.weight * exp(arma::dot(wrk, omega));
      if(grad){
        term_v.subvec(o_idx_start, o_idx_end) -= g * wrk;
        if(has_delta)
          term_v.subvec(g_idx_start, g_idx_end) -= g * Z.col(i);
      } else
        term_s -= g;
    }

    if(grad)
      term_v *= d1 * exp(fixed_effect);
    else
      term_s *= d1 * exp(fixed_effect);

    if(Y[i] > 0){
      basis(wrk, log(tstop[i]));
      if(grad){
        term_v.subvec(o_idx_start, o_idx_end) += wrk;
        if(has_delta)
          term_v.subvec(g_idx_start, g_idx_end) += Z.col(i);
      } else
        term_s += fixed_effect + arma::dot(wrk, omega);
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
    omega: coefficients for the baselin.
    omega: coefficients for the fixed effects.
    Z: design matrix.
    offsets: offsets.
    n_nodes: integer with number of Gauss-Legendre quadrature nodes.
    bound_knots: boundary knots.
    inter_knots: interior knots.
    grad: logical for whether to compute the gradient og the log-likelihood.
 */

// [[Rcpp::export(rng = false)]]
arma::vec joint_start_ll
  (arma::vec const &Y, arma::vec const &tstart, arma::vec const &tstop,
   arma::vec const &omega, arma::vec const &delta, arma::mat const &Z,
   unsigned const n_nodes, arma::vec const &bound_knots,
   arma::vec const &inter_knots, bool const grad){
  if(grad)
    return(joint_start_ll_inner<true>
             (Y, tstart, tstop, omega, delta, Z, n_nodes, bound_knots,
              inter_knots));

  return(joint_start_ll_inner<false>
           (Y, tstart, tstop, omega, delta, Z, n_nodes, bound_knots,
            inter_knots));
}
