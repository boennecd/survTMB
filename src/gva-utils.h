#ifndef GVA_UTILS_H
#define GVA_UTILS_H

#include "gaus-hermite.h"
#include "pnorm-log.h"

namespace GaussHermite {
namespace GVA {
/* TODO: make adaptive versions */

/*
 Makes an approximation of
 \begin{align*}
 l(\mu,\sigma) &=
 \frac 1{\sigma\sqrt{2\pi}}\int
 \exp \left(-\frac{(x-\mu)^2}{2\sigma^2} \right)
 \log\left(1 + \exp(x) \right)dx \\
 &\overset{x = \mu + \sqrt 2\sigma z}{=}
 \frac 1{\sqrt\pi}\int\exp(-z^2)
 \log\left(1 + \exp(\mu + \sqrt 2 \sigma z) \right)dz
 \end{align*}

 with (non-adaptive) Gauss–Hermite quadrature
*/
template<class Type>
Type mlogit_integral
  (Type const mu, Type const sigma, HermiteData<Type> const &hd){
  auto const &x = hd.x,
             &w = hd.w;
  unsigned const n_nodes = x.size();

  Type const mult_sum(sqrt(M_1_PI)),
                 mult(Type(M_SQRT2) * sigma),
            too_large(30.),
                  one(1.);
  Type out(0.);
  for(unsigned i = 0; i < n_nodes; ++i){
    Type const eta = mu + mult * x[i];
    out +=
      w[i] * CppAD::CondExpGe(eta, too_large, eta, log(one + exp(eta)));
  }

  return mult_sum * out;
}

template<class Type>
Type mlogit_integral
  (Type const mu, Type const sigma, Type const log_k,
   HermiteData<Type> const &hd){
  Type const mu_use = mu + log_k,
            sig_use = sigma;
  return mlogit_integral(mu_use, sig_use, hd);
}

/* Makes an approximation of
 l(\mu,\sigma) =
 \int\phi(x;\mu,\sigma^2)
 (-\log \Phi(x))dx

 with (non-adaptive) Gauss–Hermite quadrature
*/
template<class Type>
Type probit_integral
  (Type const mu, Type const sigma, HermiteData<Type> const &hd){
  auto const &x = hd.x,
             &w = hd.w;

  Type out(0.);
  Type const mult_sum(sqrt(M_1_PI)),
                 mult(Type(M_SQRT2) * sigma);
  for(unsigned i = 0; i < x.size(); ++i)
    out += w[i] * (-pnorm_log(mu + mult * x[i]));

  return mult_sum * out;
}

template<class Type>
Type probit_integral
  (Type const mu, Type const sigma, Type const k,
   HermiteData<Type> const &hd){
  return probit_integral(k - mu, sigma, hd);
}

} // namespace GVA
} // namespace GaussHermite

#endif
