#ifndef ORTH_POLY_H
#define ORTH_POLY_H
#include <RcppArmadillo.h>
#include <stdexcept> // invalid_argument

namespace poly {
using namespace arma;

struct orth_poly {
  vec const alpha,
            norm2,
       sqrt_norm2 = arma::sqrt(norm2);

  orth_poly(vec const &alpha, vec const &norm2):
    alpha(alpha), norm2(norm2) {
    for(size_t i = 0; i < norm2.size(); ++i)
      if(norm2[i] <= 0.)
        throw std::invalid_argument("invalid norm2");
    if(alpha.n_elem + 2L != norm2.n_elem)
      throw std::invalid_argument("invalid alpha");
  }

  /**
   behaves like predict(<poly object>, newdata) except the scale is
   different. See ::get_poly_basis().
   */
  void operator()(vec&, double const) const;
  vec operator()(double const x) const {
    vec out(get_n_basis());
    operator()(out, x);
    return out;
  };

  /**
   behaves like poly(x, degree) though the output is not scaled to have unit
   norm buth rather a norm that scales like the number of samples and there
   is an intercept. I.e. similar to

      x <- rnorm(10)
      B <- poly(x, degree = 4)
      B <- cbind(1, B * sqrt(NROW(B)))
      B
      crossprod(B)

      # same as
      B <- poly(x, degree = 4)
      attr(B, "coefs")$norm2 <- attr(B, "coefs")$norm2 / NROW(B)
      cbind(1, predict(B, x))

   The orthogonal polynomial is returned by reference.
   */
  static orth_poly get_poly_basis(vec, uword const, mat&);

  uword get_n_basis() const {
    return norm2.n_elem - 1L;
  }
};
} // namespace poly

#endif
