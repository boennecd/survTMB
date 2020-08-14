#ifndef GVA_UTILS_H
#define GVA_UTILS_H

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

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
struct mlogit_fam {
  static double const too_large;

  static double g(double const &eta) {
    return eta > too_large ? eta : log(1 + exp(eta));
  }

  template<typename T>
  static T gp(T const &eta){
    T const one(1.);

    return CppAD::CondExpGe(
      eta, T(too_large), one, one / (one + exp(-eta)));
  }
  static double gp(double const &eta) {
    return eta > too_large ? 1 : 1. / (1. + exp(-eta));
  }
};

/* atomic function to perform Gauss–Hermite quadrature.
 *
 * Args:
 *   Type: base type.
 *   Fam: class with two static function: g is the integrand and gp is the
 *        derivative.
 */
template <class Type, class Fam>
class integral_atomic : public CppAD::atomic_base<Type> {
  unsigned const n;
  HermiteData<double> const &xw_double = GaussHermiteDataCached<double>(n);
  HermiteData<Type>   const &xw_type   = GaussHermiteDataCached<Type>  (n);

public:
  integral_atomic(char const *name, unsigned const n):
  CppAD::atomic_base<Type>(name), n(n) {
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }

  /* returns a cached value to use in computations as the object must remain
   * in scope while all CppAD::ADfun functions are still in use. */
  static integral_atomic& get_cached(unsigned const);

  static double comp(
      double const mu, double const sig, HermiteData<double> const &xw){
    double out(0.);
    for(unsigned i = 0; i < xw.x.size(); ++i)
      out += xw.w[i] * Fam::g(mu + sig * xw.x[i]);
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

    ty[0] = Type(
      comp(asDouble(tx[0]), M_SQRT2 * asDouble(tx[1]), xw_double));

    /* set variable flags */
    if (vx.size() > 0) {
      bool anyvx = false;
      for (std::size_t i = 0; i < vx.size() and !anyvx; i++)
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

    Type const mu = tx[0],
             mult = Type(M_SQRT2),
              sig = mult * tx[1];
    px[0] = Type(0.);
    px[1] = Type(0.);
    for(unsigned i = 0; i < xw_type.x.size(); ++i){
      Type const node = mu + sig * xw_type.x[i],
                 term = xw_type.w[i] * Fam::gp(node);
      px[0] +=                  term;
      px[1] += mult * xw_type.x[i] * term;

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

template<class Type>
using mlogit_integral_atomic = integral_atomic<Type, mlogit_fam>;

template<class Type>
AD<Type> mlogit_integral
  (AD<Type> const mu, AD<Type> const sigma, unsigned const n){
  auto &functor = mlogit_integral_atomic<Type>::get_cached(n);

  CppAD::vector<AD<Type> > tx(2), ty(1);
  tx[0] = mu;
  tx[1] = sigma;

  functor(tx, ty);
  return ty[0];
}

inline double mlogit_integral
  (double const mu, double const sigma,  unsigned const n){
  HermiteData<double> const &xw = GaussHermiteDataCached<double>(n);
  return mlogit_integral_atomic<double>::comp(mu, M_SQRT2 * sigma, xw);
}

template<class Type>
Type mlogit_integral
  (Type const mu, Type const sigma, Type const log_k, unsigned const n){
  Type const mu_use = mu + log_k;
  return mlogit_integral(mu_use, sigma, n);
}

/* Makes an approximation of
 l(\mu,\sigma) =
 \int\phi(x;\mu,\sigma^2)
 (-\log \Phi(x))dx

 with (non-adaptive) Gauss–Hermite quadrature
*/
struct probit_fam {
  static double g(double const &eta) {
    return - atomic::Rmath::Rf_pnorm5(eta, 0, 1, 1, 1);
  }

  template<typename T>
  static T gp(T const &eta){
    T const cdf = pnorm_log(eta),
            pdf = dnorm(eta, T(0.), T(1.), 1L);
    return - exp(pdf - cdf);
  }
  static double gp(double const &eta) {
    if(eta > -10){
      auto const cdf = atomic::Rmath::Rf_pnorm5(eta, 0, 1, 1, 0);
      return - M_1_SQRT_2PI * exp(-eta * eta * .5) / cdf;
    }

    double const cdf = atomic::Rmath::Rf_pnorm5(eta, 0, 1, 1, 1),
                 pdf = - eta * eta * .5;
    return - M_1_SQRT_2PI * exp(pdf - cdf);
  }
};

template<class Type>
using probit_integral_atomic = integral_atomic<Type, probit_fam>;

template<class Type>
AD<Type> probit_integral
  (AD<Type> const mu, AD<Type> const sigma,  unsigned const n){
  auto &functor = probit_integral_atomic<Type>::get_cached(n);

  CppAD::vector<AD<Type> > tx(2), ty(1);
  tx[0] = mu;
  tx[1] = sigma;

  functor(tx, ty);
  return ty[0];
}

inline double probit_integral
  (double const mu, double const sigma, unsigned const n){
  HermiteData<double> const &xw = GaussHermiteDataCached<double>(n);
  return probit_integral_atomic<double>::comp(mu, M_SQRT2 * sigma, xw);
}

/* Notice this overload is for
 l(\mu,\sigma) =
 \int\phi(x;\mu,\sigma^2)
 (-\log \Phi(k - x))dx
 */
template<class Type>
Type probit_integral
  (Type const mu, Type const sigma, Type const k, unsigned const n){
  return probit_integral(k - mu, sigma, n);
}

} // namespace GVA
} // namespace GaussHermite

#endif
