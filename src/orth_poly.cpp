#include "orth_poly.h"

namespace poly {
void orth_poly::operator()(vec &out, double const x) const {
#ifdef DO_CHECKS
  if(out.n_elem != get_n_basis())
    throw std::invalid_argument("orth_poly::operator(): invalid out");
#endif
  out[0] = 1.;

  if(alpha.n_elem > 0L){
    out[1] = x - alpha[0];
    for(size_t c = 1; c < alpha.n_elem; c++)
      out[c + 1L] =
        (x - alpha[c]) * out[c] - norm2[c + 1L] / norm2[c] * out[c - 1L];
  }

  out /= sqrt_norm2.subvec(1L, sqrt_norm2.n_elem - 1L);
  out[0] = 1.;
}

orth_poly orth_poly::get_poly_basis(vec x, uword const degree, mat &X){
  size_t const n = x.n_elem,
              nc = degree + 1L;
  double const x_bar = mean(x);
  x -= x_bar;
  mat XX(n, nc);
  XX.col(0).ones();
  for(size_t d = 1L; d < nc; d++){
    double       * xx_new = XX.colptr(d);
    double const * xx_old = XX.colptr(d - 1);
    for(size_t i = 0; i < n; ++i, ++xx_new, ++xx_old)
      *xx_new = *xx_old * x[i];
  }

  mat R;
  if(!qr_econ(X, R, XX))
    /* TODO: can be done smarter by calling LAPACK or LINPACK directly */
    throw std::runtime_error(
        "orth_poly::get_poly_basis(): QR decomposition failed");

  for(size_t c = 0; c < nc; ++c)
    X.col(c) *= R.at(c, c);

  vec norm2(nc + 1L),
      alpha(nc - 1L);
  norm2[0] = 1.;
  for(size_t c = 0; c < nc; ++c){
    double z_sq(0),
    x_z_sq(0);
    double const *X_i = X.colptr(c);
    for(size_t i = 0; i < n; ++i, ++X_i){
      double const z_sq_i = *X_i * *X_i;
      z_sq += z_sq_i;
      if(c < degree)
        x_z_sq += x[i] * z_sq_i;
    }
    norm2[c + 1] = z_sq;
    if(c < degree)
      alpha[c] = x_z_sq / z_sq + x_bar;
  }

  orth_poly out(alpha, norm2);
  X.each_row() /= out.sqrt_norm2.subvec(1L, out.sqrt_norm2.n_elem - 1L).t();
  X.col(0).ones();
  return out;
}
} // namespace poly
