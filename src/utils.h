#ifndef SURVTMB_UTILS_H
#define SURVTMB_UTILS_H

#include "tmb_includes.h"

namespace survTMB {

inline unsigned get_rng_dim(unsigned const n_vcov_params){
  double const n(n_vcov_params);
  return std::lround(std::sqrt(.25 + 2 * n) - .5);
}

template<class Type>
unsigned get_rng_dim(vector<Type> x){
  return get_rng_dim(x.size());
}

/* Returns a covariance matix given a vector containing log standard
 * deviations, s, and lower triangular matrix elements theta. E.g. in 3x3

 \begin{align*}
 (\log \vec\sigma)^\top &= (s_1, s_2, s_3)
 & L &= \begin{pmatrix}
 1 & 0 & 0 \       \
 \theta_1 & 1 & 0 \\
 \theta_2 & \theta_3 & 1
 \end{pmatrix} \\
 \Sigma &= \text{diag}(\sigma)LL^\top\text{diag}(\sigma)
 \end{align*}
*/
template<class Type>
matrix<Type>
get_vcov_from_trian(Type const *vals, unsigned const dim){
  matrix<Type> out(dim, dim);
  out.setZero();

  /* fill the diagonal */
  for(unsigned i = 0; i < dim; ++i)
    out(i, i) = exp(*(vals + i));

  /* fill the off-diagonal */
  Type const * t = vals + dim;
  for(unsigned cl = 0; cl < dim; cl++)
    for(unsigned rw = cl + 1L; rw < dim; rw++)
      out(rw, cl) = out(rw, rw) * *t++;

  return out * out.transpose();
}

template<class Type>
matrix<Type>
get_vcov_from_trian(vector<Type> const &theta){
  unsigned const dim = get_rng_dim(theta);
  return get_vcov_from_trian(&theta[0], dim);
}
} // namespace survTMB

#endif
