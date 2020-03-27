#include "testthat-wrap.h"
#include "gva-utils.h"
#include <vector>

using namespace GaussHermite;
using namespace GaussHermite::GVA;

context("gav-utils unit tests") {
  test_that("mlogit_integral gives the correct result") {
    /*
     integrand <- function(x, mu, sigma, k, s, use_log = FALSE){
     f1  <- dnorm(x, mean = mu, sd = sigma, log = use_log)
     eta <- log(k) + s * x
     f2  <- ifelse(eta > 30, eta, log(1 + exp(eta)))
     f1 * f2
     }
     library(fastGHQuad)
     psi <- function(mu, sigma, k, s, rule){
     x      <- mu + sqrt(2) * sigma * rule$x
     f <- integrand(x, mu, sigma, k, s) * exp((x - mu)^2 / 2 / sigma^2)
     sigma * sqrt(2) * drop(rule$w %*% f)
     }
     wx <- gaussHermiteData(15L)
     dput(mapply(psi,
     mu    = c(-1,   1,  0,  2),
     sigma = c(.01,  1, 10,  5),
     k     = c(  2, .5,  5,  6),
     s     = c(  1,  1, 1, 1),
     MoreArgs = list(rule = wx)))
     */
    constexpr unsigned const n_nodes(15L);
    auto xw = GaussHermiteData(n_nodes);
    std::vector<double> x(n_nodes), w(n_nodes);
    std::vector<double> const mu = { -1,   1,  0,  2 },
                           sigma = { .01,  1, 10,  5 },
                               k = {   2, .5,  5,  6 },
                          ex_res = { 0.551456924101031, 0.969190193406026, 4.88998305310308, 4.53173698285082 };

    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = mlogit_integral(
        mu[i], sigma[i], std::log(k[i]), xw);

      expect_equal(ex_res[i], intval);
    }
  }

  test_that("probit_integral gives the correct result"){
    /*
     integrand <- function(x, mu, sigma, rho, k, s, use_log = FALSE){
    f1 <- dnorm(x, mean = mu, sd = sigma, log = use_log)
    f2 <- -pnorm(k + s * x, log.p = TRUE)
    if(use_log) f1 + log(f2) else f1 * f2
    }
    psi <- function(mu, sigma, rho, k, s)
    integrate(
      integrand, mu = mu, sigma = sigma, rho = rho, k = k, s = s,
    lower = -Inf, upper = Inf, rel.tol = 1e-10, abs.tol = 0)$value
    # the approximation should work well for these values
    dput(mapply(
        psi,
    mu    =  c(0, -1, 1, 0, -1),
    sigma =  c(1, .5, 1, 2, .7),
    k     =  c(1, -1, 2, 1, .8),
    s     = -c(1, 1, 1, 1, 1)))
    */
    constexpr unsigned const n_nodes = 40L;
    auto xw = GaussHermiteData(n_nodes);

    std::vector<double> const
        mu     = { 0, -1, 1, 0, -1 },
        sigma  = { 1, .5, 1, 2, .7 },
        k      = { 1, -1, 2, 1, .8 },
        ex_res = { 0.35934760083045, 0.771877623464394, 0.35934760083045,
                   0.942178782481152, 0.078745703164216 };
    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = probit_integral(
        mu[i], sigma[i], k[i], xw);
      expect_equal(ex_res[i], intval);
    }
  }
}
