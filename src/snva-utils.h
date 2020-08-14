#ifndef SNVA_UTILS_H
#define SNVA_UTILS_H

#include "gaus-hermite.h"
#include "pnorm-log.h"
#include "utils.h"
#include "gamma-to-nu.h"

namespace atomic {
namespace Rmath {
extern "C" {
  double Rf_dnorm4(double, double, double, int);
}
} // namespace Rmath
} // namespace atomic

namespace GaussHermite {
namespace SNVA {

/*
 Compute psi function w/ adaptive Gauss–Hermite quadrature. That is,

 \begin{align*}
 f(\sigma^2) &= \int 2 \phi(z; \sigma^2)\Phi(z)\log \Phi(z)dz \\
 &= \frac 2{\sqrt{2\pi\sigma^2}}
 \int \exp\left(-\frac {z^2}{2\sigma^2}\right)\Phi(z)\log \Phi(z)dz \\
 &= \frac 2{\sqrt{2\pi\sigma^2}}
 \int \exp\left(-\frac {z^2}
 {2\gamma^2\sigma^2/(\gamma^2 + \sigma^2)}\right)
 \underbrace{\exp\left(\frac {z^2}{2\gamma^2}\right)
 \Phi(z)\log \Phi(z)}_{f(z;\gamma)}dz \\
 &\overset{z = s \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2}}{=}
 \frac 2{\sigma\sqrt{2\pi}}
 \frac {\sqrt 2\sigma\gamma}{\sqrt{\gamma^2 + \sigma^2}}
 \int \exp\left(-s^2\right)
 f(s \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2};\gamma)ds \\
 &\approx
 \frac {2\gamma}{\sqrt{\pi(\gamma^2 + \sigma^2)}}
 \sum_{i = 1}^n
 w_i f(x_i \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2})
 \end{align*}

 \begin{align*}
 \frac{\partial f(\sigma^2)}{\partial \sigma^2} &=
 \int \frac{z^2-\sigma^2}{2\sigma^{4}}2 \phi(z; \sigma^2)\Phi(z)\log \Phi(z)dz \\
 &= \frac 1{\sqrt{2\pi\sigma^2}\sigma^{4}}
 \int (z^2 - \sigma^2)
 \exp\left(-\frac {z^2}{2\sigma^2}\right)\Phi(z)\log \Phi(z)dz \\
 &= \frac 1{\sqrt{2\pi\sigma^2}\sigma^{4}}
 \int (z^2 - \sigma^2) \exp\left(-\frac {z^2}
 {2\gamma^2\sigma^2/(\gamma^2 + \sigma^2)}\right)
 \underbrace{\exp\left(\frac {z^2}{2\gamma^2}\right)
 \Phi(z)\log \Phi(z)}_{f(z;\gamma)}dz \\
 &\overset{z = s \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2}}{=}
 \frac 1{\sigma^{5}\sqrt{2\pi}}
 \frac {\sqrt 2\sigma\gamma}{\sqrt{\gamma^2 + \sigma^2}}
 \int
 \sigma^2\left(\frac{s^2 2\gamma^2}{\gamma^2 + \sigma^2} - 1\right)
 \exp\left(-s^2\right)
 f(s \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2};\gamma)ds \\
 &\approx
 \frac 1{\sigma^2\sqrt{\pi}}
 \frac {\gamma}{\sqrt{\gamma^2 + \sigma^2}}
 \sum_{i = 1}^n
 w_i \left(\frac{x_i^2 2\gamma^2}{\gamma^2 + \sigma^2} - 1\right)
 f(x_i \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2})
 \end{align*}

 with \gamma(\sigma) = 1.
*/
template<class Type>
class entropy_term_integral : public CppAD::atomic_base<Type> {
  unsigned const n;
  HermiteData<double> const &xw_double = GaussHermiteDataCached<double>(n);
  HermiteData<Type>   const &xw_type   = GaussHermiteDataCached<Type>  (n);

  Type const zero = Type(0.),
              one = Type(1.),
              two = Type(2.),
            small = Type(std::numeric_limits<double>::epsilon() *
              std::numeric_limits<double>::epsilon()),
        type_M_PI = Type(M_PI);

public:
  entropy_term_integral(char const *name, unsigned const n):
  CppAD::atomic_base<Type>(name), n(n) {
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }

  /* returns a cached value to use in computations as the object must remain
   * in scope while all CppAD::ADfun functions are still in use. */
  static entropy_term_integral& get_cached(unsigned const);

  static double comp(double const sigma_sq, HermiteData<double> const &hd){
    double out(0.);
    double const mult_sum(M_2_SQRTPI / sqrt(sigma_sq + 1.)),
                     mult(mult_sum * sqrt(sigma_sq / M_2_PI));

    for(unsigned i = 0; i < hd.x.size(); ++i){
      double const xi = hd.x[i] * mult,
             pnrm_log = pnorm_log(xi, 0., 1.);
      out += hd.w[i] * exp(xi * xi / 2.) * exp(pnrm_log) * pnrm_log;
    }

    return mult_sum * out;
  }

  virtual bool forward(std::size_t p, std::size_t q,
                       const CppAD::vector<bool> &vx,
                       CppAD::vector<bool> &vy,
                       const CppAD::vector<Type> &tx,
                       CppAD::vector<Type> &ty){
    if(q > 0)
      return false;

    ty[0] = Type(comp(asDouble(tx[0]), xw_double));

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

    Type const  sigma_sq = tx[0] + small,
             sigma_sq_p1 = sigma_sq + one,
                mult_sum =
                  one / sigma_sq / sqrt(type_M_PI * sigma_sq_p1),
                    mult = sqrt(two * sigma_sq / sigma_sq_p1),
            mult_sq_term = two / sigma_sq_p1;

    px[0] = Type(0.);
    for(unsigned i = 0; i < n; ++i){
      Type const &x = xw_type.x[i],
                 xi = x * mult,
           pnrm_log = pnorm_log(xi, zero, one),
          integrand = exp(xi * xi / two) * exp(pnrm_log) * pnrm_log;

      px[0] += xw_type.w[i] * (x * x * mult_sq_term - one) * integrand;
    }

    px[0] *= mult_sum * py[0];

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
AD<Type> entropy_term
  (AD<Type> const sigma_sq, unsigned const n_nodes){
  auto &functor = entropy_term_integral<Type>::get_cached(n_nodes);

  CppAD::vector<AD<Type> > tx(1), ty(1);
  tx[0] = sigma_sq;

  functor(tx, ty);
  return ty[0];
}

inline double entropy_term
  (double const sigma_sq, unsigned const n_nodes){
  HermiteData<double> const &hd =
    GaussHermiteDataCached<double>(n_nodes);
  return entropy_term_integral<double>::comp(sigma_sq, hd);
}

/* Computes the approximate mode of the skew-normal distribution and the
 * Hessian at the approximation */
template<class Type>
struct mode_n_Hess {
  Type mode, Hess;
};

template<class Type>
mode_n_Hess<Type> get_SNVA_mode_n_Hess
  (Type const mu, Type const sigma, Type const rho){
  Type const one(1.),
             two(2.),
            zero(0.),
             eps(std::numeric_limits<double>::epsilon()),
           alpha = sigma * rho,
          a_sign = CppAD::CondExpLe(alpha, zero, -one, one),
              nu = Type(sqrt(M_2_PI)) * alpha / sqrt(one + alpha * alpha),
           nu_sq = nu * nu,
           gamma = Type((4. - M_PI) / 2.) * nu_sq * nu /
             pow(one - nu_sq, Type(3./2.)),
            mode = mu + sigma * (
              nu - gamma * sqrt(one - nu_sq) / two -
                a_sign / two * exp(- Type(2. * M_PI) / (
                    a_sign * alpha + eps))),
               z = rho * (mode - mu),
             phi = dnorm(z, zero, one),
             Phi = pnorm(z),
            Hess = - one / sigma / sigma - rho * rho * phi * (
              z * Phi + phi) / (Phi * Phi + eps);

  return { mode, Hess };
}

template<> inline
mode_n_Hess<double> get_SNVA_mode_n_Hess
  (double const mu, double const sigma, double const rho){
  constexpr double const one(1),
                         two(2),
                        zero(0),
                         eps(std::numeric_limits<double>::epsilon());

  double const alpha = sigma * rho,
               a_sign = alpha < 0 ? -one : one,
                   nu = sqrt(M_2_PI) * alpha / sqrt(one + alpha * alpha),
                nu_sq = nu * nu,
                gamma = (4. - M_PI) / 2. * nu_sq * nu /
                  pow(one - nu_sq, 1.5),
                 mode = mu + sigma * (
                   nu - gamma * sqrt(one - nu_sq) / two -
                     a_sign / two * exp(- 2. * M_PI / (
                         a_sign * alpha + eps))),
                    z = rho * (mode - mu),
                  phi = atomic::Rmath::Rf_dnorm4(z, zero, one, 0L),
                  Phi = atomic::Rmath::Rf_pnorm5(z, zero, one, 1L, 0L),
                 Hess = - one / sigma / sigma - rho * rho * phi * (
                   z * Phi + phi) / (Phi * Phi + eps);

  return { mode, Hess };
}

/* Makes an approximation of
 l(\mu,\sigma, \rho) =\frac{2}{\sigma\sqrt{2\pi}}
 \int \exp\left(-\frac{(z - \mu)^2}{2\sigma^2}\right)
 \Phi(\rho(z - \mu))
 \log(1 + \exp(z)) dz

 and the mode of the random effect density. The approximation seems to
 poor if there is a large skew.
 */
struct mlogit_fam {
  static double const too_large;

  template<typename T>
  static T g(T const &eta) {
    return CppAD::CondExpGe(
      eta, T(too_large), eta, log(T(1) + exp(eta)));
  }
  static double g(double const &eta) {
    return eta > too_large ? eta : log(1 + exp(eta));
  }
};

/* atomic function to perform Gauss–Hermite quadrature.
 *
 * Args:
 *   Type: base type.
 *   Fam: class with one static function g which is the integrand.
 */
template <class Type, class Fam>
class integral_atomic : public CppAD::atomic_base<Type> {
  unsigned const n;
  HermiteData<double> const &xw_double = GaussHermiteDataCached<double>(n);
  HermiteData<Type>   const &xw_type   = GaussHermiteDataCached<Type>  (n);

  Type const type_M_SQRT2 = Type(M_SQRT2),
                      one = Type(1.),
          type_M_2_SQRTPI = Type(M_2_SQRTPI);

public:
  integral_atomic(char const *name, unsigned const n):
  CppAD::atomic_base<Type>(name), n(n) {
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }

  /* returns a cached value to use in computations as the object must remain
   * in scope while all CppAD::ADfun functions are still in use. */
  static integral_atomic& get_cached(unsigned const);

  static double comp(
      double const mu, double const sig, double const rho,
      HermiteData<double> const &xw){
    auto const dvals = get_SNVA_mode_n_Hess(mu, sig, rho);
    double const xi = dvals.mode,
             lambda = 1. / sqrt(-dvals.Hess),
               mult = M_SQRT2 * lambda;

    double out(0.);
    for(unsigned i = 0; i < xw.x.size(); ++i){
      const double xx = xw.x[i],
                   zz = xi + mult * xx,
                  dif = zz - mu;
      out += xw.w[i] * Fam::g(zz) *
        exp(xx * xx - dif * dif / 2. / sig / sig) * pnorm(rho * dif);
    }
    out *= M_2_SQRTPI * lambda / sig;

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
      comp(asDouble(tx[0]), asDouble(tx[1]), asDouble(tx[2]), xw_double));

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
              sig = tx[1],
              rho = tx[2];

    auto const dvals = get_SNVA_mode_n_Hess(mu, sig, rho);
    Type const xi = dvals.mode,
           lambda = one / sqrt(-dvals.Hess),
             mult = type_M_SQRT2 * lambda;

    px[0] = Type(0.);
    px[1] = Type(0.);
    px[2] = Type(0.);
    for(unsigned i = 0; i < xw_type.x.size(); ++i){
      Type const xx = xw_type.x[i],
                 zz = xi + mult * xx,
                dif = zz - mu,
            dif_std = dif / sig,
          constants = xw_type.w[i] * Fam::g(zz) * exp(xx * xx),
               dnrm = exp(- dif_std * dif_std / 2),
              ddnrm = -dnrm,
               pnrm =         pnorm (rho * dif),
              dpnrm = atomic::dnorm1(rho * dif);
      px[0] -= constants * (
        dif_std / sig * ddnrm * pnrm + rho * dnrm * dpnrm);
      px[1] -= constants * pnrm * (dif_std * dif_std - one) / sig * ddnrm;
      px[2] += constants * dnrm * dpnrm * dif;
    }

    Type const fac =  py[0] * type_M_2_SQRTPI * lambda / sig;
    px[0] *= fac;
    px[1] *= fac;
    px[2] *= fac;

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
  (AD<Type> const mu, AD<Type> const sigma, AD<Type> const rho,
   unsigned const n_nodes){
  auto &functor = mlogit_integral_atomic<Type>::get_cached(n_nodes);

  CppAD::vector<AD<Type> > tx(3), ty(1);
  tx[0] = mu;
  tx[1] = sigma;
  tx[2] = rho;

  functor(tx, ty);
  return ty[0];
}

inline double mlogit_integral
  (double const mu, double const sigma, double const rho,
   unsigned const n_nodes){
  HermiteData<double> const &hd =
    GaussHermiteDataCached<double>(n_nodes);
  return mlogit_integral_atomic<double>::comp(mu, sigma, rho, hd);
}

template<class Type>
Type mlogit_integral
  (Type const mu, Type const sigma, Type const rho, Type const log_k,
   unsigned const n_nodes){
  Type const mu_use = mu + log_k;
  return mlogit_integral(mu_use, sigma, rho, n_nodes);
}

/* Makes an approximation of
 l(\mu,\sigma, \rho) =
 2\int\phi(z;\mu,\sigma^2)
 \Phi(\rho(z - \mu))
 (-\log \Phi(z)) dz

 and the mode of the random effect density. The approximation seems to
 poor if there is a large skew.
 */
struct probit_fam {
  template<class Type>
  static Type g(Type const &eta) {
    return - pnorm_log(eta);
  }
  static double g(double const &eta) {
    return - atomic::Rmath::Rf_pnorm5(eta, 0, 1, 1, 1);
  }
};

template<class Type>
using probit_integral_atomic = integral_atomic<Type, probit_fam>;

template<class Type>
AD<Type> probit_integral
  (AD<Type> const mu, AD<Type> const sigma, AD<Type> const rho,
   unsigned const n_nodes){
  auto &functor = probit_integral_atomic<Type>::get_cached(n_nodes);

  CppAD::vector<AD<Type> > tx(3), ty(1);
  tx[0] = mu;
  tx[1] = sigma;
  tx[2] = rho;

  functor(tx, ty);
  return ty[0];
}

inline double probit_integral
  (double const mu, double const sigma, double const rho,
   unsigned const n_nodes){
  HermiteData<double> const &hd =
    GaussHermiteDataCached<double>(n_nodes);
  return probit_integral_atomic<double>::comp(mu, sigma, rho, hd);
}

template<class Type>
Type probit_integral
  (Type const mu, Type const sigma, Type const rho, Type const k,
   unsigned const n_nodes){
  return probit_integral(k - mu, sigma, -rho, n_nodes);
}

/* The following functions maps from the input vector to the
 * parameterization used in
 * > Ormerod, J. T. (2011). Skew-normal variational approximations for
 * > Bayesian inference. Unpublished article. */
template<class Type>
struct SNVA_MD_input {
  std::vector<vector<Type> > va_mus,
                             va_rhos;
  std::vector<matrix<Type> > va_lambdas;
};

template<class Type>
class get_gamma{
  Type const c1 = Type(0.99527),
             c2 = Type(2.) * c1,
            one = Type(1.);

public:
  /* inverse of 2 * c1 * logit(gamma) - c1 */
  Type operator()(Type const gtrans) const {
    return c2 / (one + exp(-gtrans)) - c1;
  }
};

/* maps from direct parameters vector to direct parameters. See
 * get_vcov_from_trian. */

template<class Type>
SNVA_MD_input<Type> SNVA_MD_theta_DP_to_DP
  (Type const * const theta_VA, size_t const n_params,
   size_t const rng_dim){
  using survTMB::get_vcov_from_trian;
  using vecT = vector<Type>;
  using std::move;
  unsigned const n_groups = n_params / (
    rng_dim * 2L + (rng_dim * (rng_dim + 1L)) / 2L);

  SNVA_MD_input<Type> out;
  std::vector<vecT > &va_mus = out.va_mus,
                     &va_rhos = out.va_rhos;
  std::vector<matrix<Type> > &va_lambdas = out.va_lambdas;

  va_mus    .reserve(n_groups);
  va_rhos   .reserve(n_groups);
  va_lambdas.reserve(n_groups);

  Type const *t = theta_VA;
  for(unsigned g = 0; g < n_groups; ++g){
    /* insert new mu vector */
    vecT mu_vec(rng_dim);
    for(unsigned i = 0; i < rng_dim; ++i)
      mu_vec[i] = *t++;
    va_mus.emplace_back(move(mu_vec));

    /* insert new lambda matrix */
    va_lambdas.emplace_back(get_vcov_from_trian(t, rng_dim));
    t += (rng_dim * (rng_dim + 1L)) / 2L;

    matrix<Type> const &Lambda = va_lambdas.back();

    /* insert new rho vector */
    vecT rho_vec(rng_dim);
    for(unsigned i = 0; i < rng_dim; ++i)
      rho_vec[i] = *t++ / sqrt(Lambda(i, i));
    va_rhos.emplace_back(move(rho_vec));
  }

  return out;
}

/* maps from mean, Covariance matrix, and __transformed__
 * Pearson's moment coefficient of skewness to direct parameters. See
 * get_vcov_from_trian */
template<class Type>
SNVA_MD_input<Type> SNVA_MD_theta_CP_trans_to_DP
  (Type const * const theta_VA, size_t const n_params,
   size_t const rng_dim){
  using survTMB::get_vcov_from_trian;
  using vecT = vector<Type>;
  using std::move;
  unsigned const n_mu = rng_dim,
             n_lambda = (rng_dim * (rng_dim + 1L)) / 2L,
              n_per_g = n_mu + n_lambda + rng_dim,
             n_groups = n_params / n_per_g;

  SNVA_MD_input<Type> out;
  std::vector<vecT > &va_mus = out.va_mus,
                     &va_rhos = out.va_rhos;
  std::vector<matrix<Type> > &va_lambdas = out.va_lambdas;

  va_mus    .reserve(n_groups);
  va_rhos   .reserve(n_groups);
  va_lambdas.reserve(n_groups);

  get_gamma<Type> trans_g;
  Type const *t = theta_VA,
            one(1.),
      sqrt_2_pi(sqrt(2. / M_PI));

  /* intermediaries */
  vecT psi(rng_dim),
     gamma(rng_dim);

  for(unsigned g = 0; g < n_groups; ++g, t += n_per_g){
    /* get gamma parameters */
    {
      Type const *gi = t + n_mu + n_lambda;
      for(unsigned i = 0; i < rng_dim; ++i)
        gamma[i] = trans_g(*gi++);
    }

    /* Compute intermediaries and rho */
    matrix<Type> Sigma = get_vcov_from_trian(t + n_mu, rng_dim);
    vecT k(rng_dim);
    for(unsigned i = 0; i < rng_dim; ++i){
      Type const nu = gamma_to_nu(gamma[i]);
      psi[i] = sqrt(Sigma(i, i) / (one - nu * nu));
      k  [i] = nu * psi[i] / sqrt_2_pi;
    }

    /* assign mu and Lambda */
    vecT const k_sqrt_2_pi = k * sqrt_2_pi;
    vecT mu(rng_dim);
    {
      Type const *mi = t;
      for(unsigned i = 0; i < rng_dim; ++i, ++mi)
        mu[i] = *mi;
      mu -= k_sqrt_2_pi;
    }
    va_mus.emplace_back(move(mu));

    auto const k_sqrt_2_pi_mat = k_sqrt_2_pi.matrix();
    Sigma += k_sqrt_2_pi_mat * k_sqrt_2_pi_mat.transpose();

    vecT tmp = atomic::matinv(Sigma) * k;
    Type const denom = sqrt(one - vec_dot(k, tmp));

    tmp /= denom;
    va_rhos.emplace_back(move(tmp));
    va_lambdas.emplace_back(move(Sigma));
  }

  return out;
}

/* useful class to perform computations for a given
   individual */
template<class Type>
struct SNVA_cond_dens_dat {
  /* input dependend objects */
  Type const eps,
         eps_log = log(eps),
           kappa;
  unsigned const n_nodes;

  /* potentially needed constants */
  Type const mlog_2_pi_half = Type(-log(2 * M_PI) / 2.),
                        two = Type(2.);

  SNVA_cond_dens_dat
  (Type const eps, Type const kappa, unsigned const n_nodes):
  eps(eps), kappa(kappa), n_nodes(n_nodes) { }
};

#define SNVA_COND_DENS_ARGS                                      \
  Type const &eta_fix, Type const &etaD_fix,                     \
  Type const &event, Type const &va_mu, Type const &va_sd,       \
  Type const &va_rho, Type const &va_d, Type const &va_var,      \
  Type const &dist_mean, Type const &dist_var

/* computes the conditional density term for the PH (log-log) link
 * function */
template<class Type>
struct ph final : public SNVA_cond_dens_dat<Type> {
  using SNVA_cond_dens_dat<Type>::SNVA_cond_dens_dat;

  Type operator()(SNVA_COND_DENS_ARGS) const {
    Type const h = etaD_fix * exp(eta_fix + dist_mean),
               H = this->two * exp(
                 eta_fix + va_mu + va_var / this->two) * pnorm(va_d),
          if_low = event * this->eps_log - H - h * h * this->kappa,
          if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, this->eps, if_ok, if_low);
  }
};

/* computes the conditional density term for the PO (-logit) link
 * function */
template<class Type>
struct po final : public SNVA_cond_dens_dat<Type> {
  using SNVA_cond_dens_dat<Type>::SNVA_cond_dens_dat;

  Type operator()(SNVA_COND_DENS_ARGS) const {
    Type const H = mlogit_integral(
      va_mu, va_sd, va_rho, eta_fix, this->n_nodes),
               h = etaD_fix * exp(eta_fix + dist_mean - H),
          if_low = event * this->eps_log - H - h * h * this->kappa,
          if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, this->eps, if_ok, if_low);
  }
};

/* computes the conditional density term for the probit (-probit) link
 * function */
template<class Type>
struct probit final : public SNVA_cond_dens_dat<Type> {
  using SNVA_cond_dens_dat<Type>::SNVA_cond_dens_dat;

  Type operator()(SNVA_COND_DENS_ARGS) const {
    Type const H = probit_integral(
        va_mu, va_sd, va_rho, -eta_fix, this->n_nodes),
            diff = (eta_fix + dist_mean),
               h = etaD_fix * exp(
                 this->mlog_2_pi_half - diff * diff / this->two -
                   dist_var / this->two + H),
            if_low = event * this->eps_log - H - h * h * this->kappa,
            if_ok  = event * log(h)  - H;

    return CppAD::CondExpGe(h, this->eps, if_ok, if_low);
  }
};

#undef SNVA_COND_DENS_ARGS

} // namespace SNVA
} // namespace GaussHermite

#endif
