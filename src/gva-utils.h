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

template <class Type>
class mlogit_integral_atomic : public CppAD::atomic_base<Type> {
  static constexpr unsigned const n = 20L; /* TODO: handle variable n */
  static constexpr double const too_large = 30.;

  template<typename T>
  static T g(T const &eta) {
    return CppAD::CondExpGe(
        eta, T(too_large), eta, log(T(1) + exp(eta)));
  }

  template<typename T>
  static T gp(T const &eta){
    T const one(1.);

    return CppAD::CondExpGe(
      eta, T(too_large), one, one / (one + exp(-eta)));
  }

public:
  mlogit_integral_atomic(char const *name):
  CppAD::atomic_base<Type>(name) {
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }

  static double comp(
      double const mu, double const sig, HermiteData<double> const &xw){
    double out(0.);
    for(unsigned i = 0; i < xw.x.size(); ++i)
      out += xw.w[i] * g(mu + sig * xw.x[i]);
    out *= sqrt(M_1_PI);

    return out;
  }

  virtual bool forward(std::size_t p, std::size_t q,
                       const CppAD::vector<bool> &vx,
                       CppAD::vector<bool> &vy,
                       const CppAD::vector<Type> &tx,
                       CppAD::vector<Type> &ty){
    if(q > 0)
      return false;

    HermiteData<double> const &xw = GaussHermiteDataCached(n);
    ty[0] = Type(
      comp(asDouble(tx[0]), M_SQRT2 * asDouble(tx[1]), xw));

    /* set variable flags */
    if (vx.size() > 0) {
      bool anyvx = false;
      for (std::size_t i = 0; i < vx.size(); i++)
        anyvx |= vx[i];
      for (std::size_t i = 0; i < vy.size(); i++)
        vy[i] = anyvx;
    }

    return true;
  }

  virtual bool reverse(std::size_t q, const CppAD::vector<Type> &tx,
                       const CppAD::vector<Type> &ty,
                       CppAD::vector<Type> &px,
                       const CppAD::vector<Type> &py){
    if(q > 0)
      return false;

    /* TODO: switch to <Type> */
    HermiteData<double> const &xw = GaussHermiteDataCached(n);
    Type const mu = tx[0],
             mult = Type(M_SQRT2),
              sig = mult * tx[1];
    px[0] = Type(0.);
    px[1] = Type(0.);
    for(unsigned i = 0; i < xw.x.size(); ++i){
      Type const xx = Type(xw.x[i]),
               node = mu + sig * xx,
               term = Type(xw.w[i]) * gp(node);
      px[0] +=             term;
      px[1] += mult * xx * term;

    }

    px[0] *= Type(sqrt(M_1_PI)) * py[0];
    px[1] *= Type(sqrt(M_1_PI)) * py[0];

    return true;
  }

  virtual bool rev_sparse_jac(size_t q, const CppAD::vector<bool>& rt,
                              CppAD::vector<bool>& st) {
    bool anyrt = false;
    for (std::size_t i = 0; i < rt.size(); i++)
      anyrt |= rt[i];
    for (std::size_t i = 0; i < st.size(); i++)
      st[i] = anyrt;
    return true;
  }
};

/* TODO: get rid of HermiteData<AD<Type> > arg */
template<class Type>
AD<Type> mlogit_integral
  (AD<Type> const mu, AD<Type> const sigma, HermiteData<AD<Type> > const &hd){
  static mlogit_integral_atomic<Type> functor("mlogit");

  CppAD::vector<AD<Type> > tx(2), ty(1);
  tx[0] = mu;
  tx[1] = sigma;

  functor(tx, ty);
  return ty[0];
}

double mlogit_integral
  (double const mu, double const sigma, HermiteData<double> const &hd);

template<class Type>
Type mlogit_integral
  (Type const mu, Type const sigma, Type const log_k,
   HermiteData<Type> const &hd){
  Type const mu_use = mu + log_k;
  return mlogit_integral(mu_use, sigma, hd);
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
