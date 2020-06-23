#ifndef GSM_H
#define GSM_H

#define INCLUDE_RCPP
#include "tmb_includes.h"
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace gsm_objs {
/** evalutes the dot product between an arma::vec and double range without
 bounds checks. */
inline double arma_dot(const arma::vec &x, double const *y){
  double out(0.);
  for(auto &xi : x)
    out += xi * *y++;
  return out;
}

/** adds a double range times a scalar to an arma::vec without bounds
 checks. */
inline void arma_inplace_add(arma::vec &x, double const f, double const *y){
  for(auto &xi : x)
    xi += f * *y++;
}

/** abstract base class to return to R. */
class gsm_base {
public:
  virtual double log_likelihood
  (arma::vec const&, arma::vec const&) const = 0;
  virtual arma::vec grad(arma::vec const&, arma::vec const&) const = 0;
  virtual arma::mat hess(arma::vec const&, arma::vec const&) const = 0;

  virtual ~gsm_base() = default;
};

/**
 template class to do the computation of the log-likelihood and the
 gradient. The template parameter needs a constructor which takes the linear
 predictor and has the following member functions:
   g_log: log of the inverse link function.
   gp: derivative of the inverse link function.
   gp_g: ratio of the derivative and the the inverse link function.
   gpp_gp: ratio of twice-derivative and the derivative of the inverse
           link function.
   gpp: twice-derivative of the inverse link function.
   d_gp_g: derivative of gp_g.
   d_gpp_gp: derivative of gpp_gp.

 This allows one to cache certian values. The design matrices needs to be
 [# coefficients] x [# observations]. Z is for the time invariant
 covariates.
 */
template<class Family>
class gsm final : public gsm_base {
  arma::mat const X, XD, Z;
  arma::vec const y;
  size_t const n = X.n_cols,
             n_b = X.n_rows,
             n_g = Z.n_rows;
  double const eps, kappa,
           eps_log = log(eps);
  unsigned const n_threads;
  arma::vec const offset_eta, offset_etaD;

  void check_params(arma::vec const &beta, arma::vec const &gamma) const {
    if(beta.n_elem != n_b)
      throw std::invalid_argument("gsm: invalid beta");
    else if(gamma.n_elem != n_g)
      throw std::invalid_argument("gsm: invalid gamma");
  }

public:
  gsm(arma::mat const &X, arma::mat const &XD, arma::mat const &Z,
      arma::vec const &y, double const eps, double const kappa,
      unsigned const n_threads, arma::vec const &offset_eta,
      arma::vec const &offset_etaD):
  X(X), XD(XD), Z(Z), y(y), eps(eps), kappa(kappa), n_threads(n_threads),
  offset_eta(offset_eta), offset_etaD(offset_etaD) {
    /* checks */
    if(XD.n_rows != n_b or XD.n_cols != n)
      throw std::invalid_argument("gsm: invalid XD");
    else if(Z.n_cols != n)
      throw std::invalid_argument("gsm: invalid Z");
    else if(y.n_elem != n)
      throw std::invalid_argument("gsm: invalid y");
    else if(kappa < 0)
      throw std::invalid_argument("gsm: invalid kappa");
    else if(eps < 0)
      throw std::invalid_argument("gsm: invalid eps");
    else if(n_threads < 1)
      throw std::invalid_argument("gsm: invalid n_threads");
    else if(offset_eta.n_elem != n)
      throw std::invalid_argument("gsm: invalid offset_eta");
    else if(offset_etaD.n_elem != n)
      throw std::invalid_argument("gsm: invalid offset_etaD");

#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#endif
  }

  double log_likelihood
  (arma::vec const &beta, arma::vec const &gamma) const {
    check_params(beta, gamma);

    double out(0.);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) schedule(static) \
    reduction(+:out)
#endif
    for(size_t i = 0; i < n; ++i){
      double const * const xi  = X .colptr(i),
                   * const xdi = XD.colptr(i),
                   * const zi  = Z .colptr(i);
      double const eta =
        arma_dot(beta, xi) + arma_dot(gamma, zi) + offset_eta[i];
      Family fam(eta);
      double const eta_p = arma_dot(beta, xdi) + offset_etaD[i],
                     haz = -fam.gp_g() * eta_p;
      bool const valid = haz > eps;


      if(y[i] > 0)
        if(__builtin_expect(valid, 1))
          out += log(-fam.gp() * eta_p);
        else
          out += eps_log + fam.g_log();
      else
        out += fam.g_log();

      if(__builtin_expect(valid, 1))
        continue;

      double const delta = haz - eps;
      out -= kappa * delta * delta;
    }

    return out;
  }

  arma::vec grad
  (arma::vec const &beta, arma::vec const &gamma) const {
    check_params(beta, gamma);
    arma::vec out(n_b + n_g, arma::fill::zeros),
               db(out.begin()      , n_b, false),
               dg(out.begin() + n_b, n_g, false);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
    {
#endif
    arma::vec db_loc(n_b, arma::fill::zeros),
              dg_loc(n_g, arma::fill::zeros);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(size_t i = 0; i < n; ++i){
      double const * const xi  = X .colptr(i),
                   * const xdi = XD.colptr(i),
                   * const zi  = Z .colptr(i);
      double const eta =
          arma_dot(beta, xi) + arma_dot(gamma, zi) + offset_eta[i];
      Family fam(eta);
      double const eta_p = arma_dot(beta, xdi) + offset_etaD[i],
                     haz = -fam.gp_g() * eta_p;
      bool const valid = haz > eps;

      if(y[i] > 0){
        if(__builtin_expect(valid, 1)){
          double const f = fam.gpp_gp();
          arma_inplace_add(db_loc, f         , xi);
          arma_inplace_add(dg_loc, f         , zi);
          arma_inplace_add(db_loc, 1. / eta_p, xdi);

        } else {
          double const fac = fam.gp_g();
          arma_inplace_add(db_loc, fac, xi);
          arma_inplace_add(dg_loc, fac, zi);

        }
      } else {
        double const fac = fam.gp_g();
        arma_inplace_add(db_loc, fac, xi);
        arma_inplace_add(dg_loc, fac, zi);

      }

      if(__builtin_expect(valid, 1))
        continue;

      double const fac = -2. * kappa * (haz - eps),
                    f1 = -fac * fam.d_gp_g() * eta_p,
                    f2 = -fac * fam.gp_g();
      arma_inplace_add(db_loc, f1, xi);
      arma_inplace_add(dg_loc, f1, zi);
      arma_inplace_add(db_loc, f2, xdi);
    }

#ifdef _OPENMP
#pragma omp critical
      {
#endif
    db += db_loc;
    dg += dg_loc;
#ifdef _OPENMP
      }
    }
#endif

    return out;
  }

  arma::mat hess
  (arma::vec const &beta, arma::vec const &gamma) const {
    check_params(beta, gamma);
    arma::mat bm(n_b, n_b, arma::fill::zeros),
              gm(n_g, n_g, arma::fill::zeros),
             gbm(n_g, n_b, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
  {
#endif
    arma::mat bm_loc(n_b, n_b, arma::fill::zeros),
              gm_loc(n_g, n_g, arma::fill::zeros),
             gbm_loc(n_g, n_b, arma::fill::zeros);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(size_t i = 0; i < n; ++i){
      double const * const xi  = X .colptr(i),
                   * const xdi = XD.colptr(i),
                   * const zi  = Z .colptr(i);
      double const eta =
        arma_dot(beta, xi) + arma_dot(gamma, zi) + offset_eta[i];
      Family fam(eta);
      double const eta_p = arma_dot(beta, xdi) + offset_etaD[i],
                     haz = -fam.gp_g() * eta_p;
      bool const valid = haz > eps;

      if(y[i] > 0){
        if(__builtin_expect(valid, 1)){
          double const f = fam.d_gpp_gp();
          /* TODO: do something smarter */
          bm_loc  += (f * X.col(i)) * X.col(i).t();
          gm_loc  += (f * Z.col(i)) * Z.col(i).t();
          gbm_loc += (f * Z.col(i)) * X.col(i).t();
          bm_loc  -= (XD.col(i) / (eta_p * eta_p)) * XD.col(i).t();

        }
      } else {
        /* TODO: do something smarter */
        double const f = fam.d_gp_g();
        bm_loc  += (f * X.col(i)) * X.col(i).t();
        gm_loc  += (f * Z.col(i)) * Z.col(i).t();
        gbm_loc += (f * Z.col(i)) * X.col(i).t();

      }
    }

#ifdef _OPENMP
#pragma omp critical
      {
#endif
    bm  += bm_loc;
    gm  += gm_loc;
    gbm += gbm_loc;
#ifdef _OPENMP
      }
    }
#endif

    arma::mat out(n_b + n_g, n_b + n_g);
    if(n_b > 0)
      out.submat(0  , 0  , n_b - 1      , n_b - 1      ) = bm;
    if(n_g > 0)
      out.submat(n_b, n_b, n_b + n_g - 1, n_b + n_g - 1) = gm;
    if(n_b > 0 and n_g > 0)
      out.submat(n_b, 0  , n_b + n_g - 1, n_b - 1      ) = gbm;
    return arma::symmatl(out);
  }
};

/** probit link function. */
struct gsm_probit {
  double const eta,
          dnrm_log = Rf_dnorm4(-eta, 0, 1, 1),
          pnrm_log = Rf_pnorm5(-eta, 0, 1, 1, 1);

  gsm_probit(double const eta): eta(eta) { }

  double g_log() const;
  double gp() const;
  double gp_g() const;
  double gpp_gp() const;
  double gpp() const;
  double d_gp_g() const;
  double d_gpp_gp() const;
};

/** PH link function. */
class gsm_ph {
  double const eta,
           exp_eta = exp(eta);

public:
  gsm_ph(double const eta): eta(eta) { }

  double g_log() const;
  double gp() const;
  double gp_g() const;
  double gpp_gp() const;
  double gpp() const;
  double d_gp_g() const;
  double d_gpp_gp() const;
};

/** negative-logit link function. */
class gsm_logit {
  double const eta,
           exp_eta = exp(eta),
        exp_eta_p1 = 1. + exp_eta;

public:
  gsm_logit(double const eta): eta(eta) { }

  double g_log() const;
  double gp() const;
  double gp_g() const;
  double gpp_gp() const;
  double gpp() const;
  double d_gp_g() const;
  double d_gpp_gp() const;
};
} // namespace gsm_objs

#endif
