#define INCLUDE_RCPP
#include "tmb_includes.h"
#include "fastgl.h"
#include "bases-wrapper.h"

template<class Basis>
arma::vec joint_start_ll_inner
  (arma::vec const &Y, arma::vec const &tstart, arma::vec const &tstop,
   arma::vec const &omega, arma::vec const &delta, arma::mat const &Z,
   unsigned const n_nodes, arma::vec const &coefs, bool const grad,
   bool const use_log){
  auto const basis = get_basis<Basis>(coefs);
  bool const has_b = static_cast<bool>(basis.get());
  size_t const dim_basis = has_b ? basis->get_n_basis() : 0L;
  size_t const n_out = grad ? dim_basis  + delta.size() : 1L;
  arma::vec out(n_out, arma::fill::zeros);
  auto const &xw = fastgl::GLPairsCached<double>(n_nodes);
  arma::vec wrk(dim_basis), term_v(n_out);

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
  else if(omega.n_elem < 1L or omega.n_elem != basis->get_n_basis())
    throw std::invalid_argument("joint_start_ll: invalid omega");
  else if(n_nodes == 0L)
    throw std::invalid_argument("joint_start_ll: invalid n_nodes");
#endif

  auto eval_basis = [&](double x){
    if(use_log)
      x = log(x);
    if(has_b)
      basis->operator()(wrk, x);
  };

  for(size_t i = 0; i < n; ++i){
    double term_s(0.);
    if(grad)
      term_v.zeros();
    double const d1 = (tstop[i] - tstart[i]) / 2.,
                 d2 = (tstop[i] + tstart[i]) / 2.,
       fixed_effect = has_delta ? arma::dot(Z.col(i), delta) : 0.;
    for(auto const &xwi : xw){
      double const node = d1 * xwi.x + d2;
      eval_basis(node);
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
      eval_basis(tstop[i]);
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
    coefs: input for the basis.
    grad: logical for whether to compute the gradient og the log-likelihood.
    use_log: logical for whether to use log(time) in the basis.
    basis_type: string with the basis type.
 */

// [[Rcpp::export(rng = false)]]
arma::vec joint_start_ll
  (arma::vec const &Y, arma::vec const &tstart, arma::vec const &tstop,
   arma::vec const &omega, arma::vec const &delta, arma::mat const &Z,
   unsigned const n_nodes, arma::vec const &coefs,
   bool const grad, bool const use_log, std::string const basis_type){
  if     (basis_type == "ns")
    return(joint_start_ll_inner<splines::ns>
             (Y, tstart, tstop, omega, delta, Z, n_nodes, coefs, grad,
              use_log));
  else if(basis_type == "poly")
    return(joint_start_ll_inner<poly::orth_poly>
             (Y, tstart, tstop, omega, delta, Z, n_nodes, coefs, grad,
              use_log));

  throw std::invalid_argument(
      "joint_start_ll: 'basis_type' not implemented");
  return arma::vec();
}
