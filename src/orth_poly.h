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

  /** behaves like predict(<poly object>, newdata). */
  void operator()(vec&, double const) const;
  vec operator()(double const x) const {
    vec out(get_n_basis());
    operator()(out, x);
    return out;
  };

  /**
   behaves like poly(x, degree). The orthogonal polynomial is returned by
   reference.
   */
  static orth_poly get_poly_basis(vec, uword const, mat&);

  uword get_n_basis() const {
    return norm2.n_elem - 1L;
  }
};
} // namespace poly

#endif
