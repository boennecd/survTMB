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

  void check_params(arma::vec const &beta, arma::vec const &gamma) const {
    if(beta.n_elem != n_b)
      throw std::invalid_argument("gsm: invalid beta");
    else if(gamma.n_elem != n_g)
      throw std::invalid_argument("gsm: invalid gamma");
  }

public:
  gsm(arma::mat const &X, arma::mat const &XD, arma::mat const &Z,
      arma::vec const &y, double const eps, double const kappa,
      unsigned const n_threads):
  X(X), XD(XD), Z(Z), y(y), eps(eps), kappa(kappa), n_threads(n_threads) {
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
      double const eta = arma_dot(beta, xi) + arma_dot(gamma, zi);
      Family fam(eta);

      if(y[i] > 0){
        double const eta_p = arma_dot(beta, xdi),
                         h = -fam.gp() * eta_p;
        if(__builtin_expect(h > eps, 1))
          /*   valid h */
          out += log(h);
        else
          /* invalid h */
          out += eps_log - h * h * kappa;

      } else
        out += fam.g_log();
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
      double const eta = arma_dot(beta, xi) + arma_dot(gamma, zi);
      Family fam(eta);

      if(y[i] > 0){
        double const eta_p = arma_dot(beta, xdi),
                         h = -fam.gp() * eta_p;

        if(__builtin_expect(h > eps, 1)){
          /*   valid h */
          double const f = fam.gpp_gp();
          arma_inplace_add(db_loc, f         , xi);
          arma_inplace_add(dg_loc, f         , zi);
          arma_inplace_add(db_loc, 1. / eta_p, xdi);

        } else {
          /* invalid h */
          double const fpp = 2. * kappa * fam.gpp(),
                        fp = 2. * kappa * fam. gp();
          arma_inplace_add(db_loc, -fpp, xi);
          arma_inplace_add(dg_loc, -fpp, zi);
          arma_inplace_add(db_loc, - fp, xdi);

        }
      } else {
        double const f = fam.gp_g();
        arma_inplace_add(db_loc, f, xi);
        arma_inplace_add(dg_loc, f, zi);

      }
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
};

/** probit link function. */
struct gsm_probit {
  double const eta;

  gsm_probit(double const eta): eta(eta) { }

  double g_log() const;
  double gp() const;
  double gp_g() const;
  double gpp_gp() const;
  double gpp() const;
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
};
} // namespace gsm_objs

#endif
