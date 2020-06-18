#ifndef BASES_WRAPPER
#define BASES_WRAPPER
#include "splines.h"
#include "orth_poly.h"
#include <memory>

/** wraps various basis given a single input vector */
template<class Basis>
std::unique_ptr<Basis> get_basis(arma::vec const&);

template<>
inline std::unique_ptr<splines::ns>
get_basis(arma::vec const &coefs){
  using splines::ns;
  if(coefs.n_elem < 2L)
    return std::unique_ptr<ns>();

  arma::vec bk(2L);
  bk[0L] = asDouble(coefs[0L]);
  bk[1L] = asDouble(coefs[coefs.n_elem - 1L]);

  arma::vec ik(coefs.size() - 2L);
  for(int i = 1L; i < coefs.n_elem - 1L; ++i)
    ik[i - 1L] = coefs[i];

  return std::unique_ptr<ns>(new ns(bk, ik, true));
}

template<>
inline std::unique_ptr<poly::orth_poly>
get_basis(arma::vec const &coefs) {
  using pol = poly::orth_poly;
  if(coefs.n_elem < 2)
    return std::unique_ptr<pol>();

  size_t n_alpha = coefs.n_elem / 2L - 1L;
  arma::vec alpha(n_alpha), norm2(coefs.n_elem - n_alpha);
  size_t i = 0;
  for(; i < n_alpha; ++i)
    alpha[i] = coefs[i];
  for(; i < coefs.n_elem; ++i)
    norm2[i - n_alpha] =  asDouble(coefs[i]);

  return std::unique_ptr<pol>(new pol(alpha, norm2));
}

#endif
