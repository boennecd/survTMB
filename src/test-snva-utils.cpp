#include "testthat-wrap.h"
#include "snva-utils.h"
#include <vector>

using namespace GaussHermite;
using namespace GaussHermite::SNVA;

context("snva-utils unit tests") {
  test_that("entropy_term gives the correct result") {
    /*
     psi <- function(sigma)
     ifelse(
     sigma < .Machine$double.eps^(1/2),
     -log(2),
     sapply(sigma, function(s)
     integrate(function(x)
     2 * dnorm(x, sd = s) * pnorm(x) * pnorm(x, log.p = TRUE),
     -Inf, Inf, rel.tol = .Machine$double.eps^(3/4))$value))
     dput(psi(0:4))
     */
    std::vector<double> sigmas = { 0., 1., 2., 100. };
    std::vector<double> expect = {
      -0.693147180559945, -0.5, -0.319885055354923, -0.00720608477031823 };
    constexpr unsigned const n_nodes = 20L;
    auto xw = GaussHermite::GaussHermiteData(n_nodes);

    for(unsigned i = 0L; i < 4L; ++i){
      double res = entropy_term(sigmas[i] * sigmas[i], xw);
      expect_equal(expect[i], res);
    }

    {
      double res = entropy_term(1e-8, xw);
      expect_equal(-std::log(2), res);
    }
    {
      double res = entropy_term(0., xw);
      expect_equal(-std::log(2), res);
    }
  }

  test_that("mlogit_integral gives the correct result") {
    /*
     integrand <- function(x, mu, sigma, rho, k, s, use_log = FALSE){
     f1 <- dnorm(x, mean = mu, sd = sigma, log = use_log)
     f2 <- 2 * pnorm(rho * (x - mu), log.p = use_log)
     eta <- log(k) + s * x
     f3  <- ifelse(eta > 30, eta, log(1 + exp(eta)))
     if(use_log) f1 + f2 + log(f3) else f1 * f2 * f3
     }
     psi <- function(mu, sigma, rho, k, s)
     integrate(
     integrand, mu = mu, sigma = sigma, rho = rho, k = k, s = s,
     lower = -Inf, upper = Inf, rel.tol = 1e-10, abs.tol = 0)$value
     # the approximation should work well for these values
     dput(mapply(
     psi,
     mu    = c(0, -1, 1, 0, -1),
     sigma = c(1, .5, 1, 2, .7),
     rho   = c(0, -1, 1, 1, -2),
     k     = c(1, .3, 2, 1, .8),
     s     = c(1, 1, 1, 1, 1)))
     */
    constexpr unsigned const n_nodes = 40L;
    auto xw = GaussHermite::GaussHermiteData(n_nodes);
    std::vector<double> const
        mu     = { 0, -1, 1, 0, -1 },
        sigma  = { 1, .5, 1, 2, .7 },
        rho    = { 0, -1, 1, 1, -2 },
        k      = { 1, .3, 2, 1, .8 },
        ex_res = { 0.80605918334744, 0.0969646948660297, 2.38729257810521, 1.78136403451249,
                   0.189741259960398 };

    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = mlogit_integral(
        mu[i], sigma[i], rho[i], std::log(k[i]), xw);
      expect_equal(ex_res[i], intval);
    }
  }

  test_that("probit_integral gives the correct result") {
    /*
     integrand <- function(x, mu, sigma, rho, k, s, use_log = FALSE){
     f1 <- dnorm(x, mean = mu, sd = sigma, log = use_log)
     f2 <- pnorm(rho * (x - mu), log.p = use_log)
     f3 <- -2 * pnorm(k + s * x, log.p = TRUE)
     f1 * f2 * f3
     }
     psi <- function(mu, sigma, rho, k, s)
     integrate(
     integrand, mu = mu, sigma = sigma, rho = rho, k = k, s = s,
     lower = -Inf, upper = Inf, rel.tol = 1e-10, abs.tol = 0)$value
     # the approximation should work well for these values
     dput(mapply(
     psi,
     mu    =  c(0, -1,  1,  0,  -1),
     sigma =  c(1, .5,  1,  2,  .7),
     rho   =  c(0, -1,  1,  1,  -2),
     k     =  c(1, .3, -2,  1, -.8),
     s     = -c(1,  1,  1,  1,   1)))
     */
    unsigned const n_nodes = 40L;
    auto xw = GaussHermite::GaussHermiteData(n_nodes);
    std::vector<double> const
        mu     = { 0, -1,  1,  0,  -1 },
        sigma  = { 1, .5,  1,  2,  .7 },
        rho    = { 0, -1,  1,  1,  -2 },
        k      = { 1, .3, -2,  1, -.8 },
        ex_res = { 0.35934760083045, 0.0981790817582183, 8.93028198953091, 1.79323321188706,
                   0.360643279313754 };
    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = probit_integral(
        mu[i], sigma[i], rho[i], k[i], xw);

      expect_equal(ex_res[i], intval);
    }
  }

  test_that("SNVA_MD_theta_DP_to_DP gives the correct result") {
    /*
     dput(mu <- c(4, 3, 8))
     dput(Sig <- matrix(c(1, 1, 2, 1, 4, 3, 2, 3, 6), 3))
     dput(rho <- c(5, 7, 9))
     ch <- t(chol(Sig))
     log_sigs <- log(diag(ch))
     L <- diag(diag(ch)^(-1)) %*% ch
     dput(c(mu, log_sigs, L[lower.tri(L)], rho))
     */
    using Type = double;
    vector<Type> theta(12);
    theta << 4, 3, 8, 0, 0.549306144334055, 0.255412811882995,
             0.577350269189626,
             1.54919333848297, 0.447213595499958, 5, 7, 9;

    auto input = SNVA_MD_theta_DP_to_DP(theta, 3L);
    std::vector<double> const mu = { 4, 3, 8 },
                             Sig = { 1, 1, 2, 1, 4, 3, 2, 3, 6 },
                             rho = { 5, 7, 9 };

    expect_true(input.va_mus.size() == 1L);

    auto do_test = [&](auto &ex, auto &re){
      expect_true(ex.size() == (unsigned)re.size());
      for(unsigned i = 0; i < ex.size(); ++i)
        expect_equal(ex[i], *(re.data() + i));
    };

    do_test(mu , input.va_mus[0]);
    do_test(Sig, input.va_lambdas[0]);
    do_test(rho, input.va_rhos[0]);
  }
  {
    /*
     cp_to_dp <- function(mu, Sigu, gamma){
     Sig <- matrix(0., length(mu), length(mu))
     Sig[upper.tri(Sig, TRUE)] <- Sigu
     Sig <- crossprod(Sig)

     gamma <- 2 * 0.99527 / (1 + exp(-gamma)) - 0.99527
     cv <- 2 * abs(gamma) / (4 - pi)
     nu <- ifelse(gamma < 0, -1, 1) * cv^(1/3) / sqrt(1 + cv^(2/3))
     omegas <- sqrt(diag(Sig) / (1 - nu^2))
     rhos <- sqrt(pi) * nu / sqrt(2 - pi * nu * nu) / omegas
     onu <- omegas * nu
     Lambda <- Sig + onu %o% onu

     list(xi = mu - onu, Lambda = Lambda, rhos = rhos)
     }
     dput(cp_to_dp(
     mu = c(-1, 0, 1),
     Sigu = Sigu <- c(4, 2, 1, 1, 2, 3),
     gamma = c(-2, .5, 0)))

     ch <- matrix(0., 3L, 3L)
     ch[lower.tri(ch, TRUE)] <- Sigu
     log_sig <- log(diag(ch))
     L <- diag(diag(ch)^(-1)) %*% ch
     dput(c(log_sig, L[lower.tri(L)]))
     */

    using Type = double;
    vector<Type> theta(12);
    theta << -1., 0., 1.,
             1.38629436111989, 0, 1.09861228866811, 2, 0.333333333333333,
             0.666666666666667,
             -2., .5, 0.;

    auto input = SNVA_MD_theta_CP_trans_to_DP(theta, 3L);
    std::vector<double> const
      mu  = { 3.83496892053337, -1.85176038792831, 1 },
      Sig = { 39.3769244625236, -0.953203923908205, 4, -0.953203923908205,
              8.42901653430041, 4, 4, 4, 14 },
      rho = { -0.592481457389101, 0.458273315832141, 0 };

    expect_true(input.va_mus.size() == 1L);

    auto do_test = [&](auto &ex, auto &re){
      expect_true(ex.size() == (unsigned)re.size());
      for(unsigned i = 0; i < ex.size(); ++i)
        expect_equal(ex[i], *(re.data() + i));
    };

    do_test(mu , input.va_mus[0]);
    do_test(Sig, input.va_lambdas[0]);
    do_test(rho, input.va_rhos[0]);
  }
}
