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
     wx <- gaussHermiteData(20L)
     dput(mapply(psi,
     mu    = c(-1,   1,  0,  2),
     sigma = c(.01,  1, 10,  5),
     k     = c(  2, .5,  5,  6),
     s     = c(  1,  1, 1, 1),
     MoreArgs = list(rule = wx)))
     dpsi <- function(mu, sigma, k){
     func <- function(x)
     psi(x[1], x[2], x[3], s = 1, rule = wx)
     numDeriv::jacobian(func, c(mu, sigma, k), method.args = list(eps = 1e-10))
     }
     dput(mapply(dpsi,
     mu    = c(-1,   1,  0,  2),
     sigma = c(.01,  1, 10,  5),
     k     = c(  2, .5,  5,  6)))
     */
    constexpr unsigned const n_nodes(20L);
    std::vector<double> x(n_nodes), w(n_nodes);
    std::vector<double> const
      mu     = { -1,   1,  0,  2 },
      sigma  = { .01,  1, 10,  5 },
      k      = {   2, .5,  5,  6 },
      ex_res = { 0.551456924101031, 0.969190193359827, 4.91705177807963, 4.53659454642133 },
      derivs = { 0.42388497395228, 0.00244200473481651, 0.211942486987524,
                 0.563103209889903, 0.203707020761764, 1.12620641974815, 0.53359761289107,
                 0.394524085551894, 0.106706563689337, 0.762029011954127, 0.293155791137677,
                 0.127004835341413 };

    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = mlogit_integral(
        mu[i], sigma[i], std::log(k[i]), n_nodes);

      expect_equal(ex_res[i], intval);
    }

    /* check Jacobian */
    using ADd = AD<double>;
    vector<ADd > a(3), b(1);
    a[0] = ADd(mu[0]);
    a[1] = ADd(sigma[0]);
    a[2] = ADd(k[0]);

    CppAD::Independent(a);
    b[0] = mlogit_integral(a[0], a[1], log(a[2]), n_nodes);
    CppAD::ADFun<double> func(a, b);

    auto d = derivs.cbegin();
    vector<double> aa(3), weight(1);
    weight[0] = 1.;
    double const eps_deriv = std::pow(
      std::numeric_limits<double>::epsilon(), 1./ 4.);
    for(unsigned i = 0; i < mu.size(); ++i){
      aa[0] = mu[i];
      aa[1] = sigma[i];
      aa[2] = k[i];

      auto yy = func.Forward(0, aa);
      expect_equal(ex_res[i], yy[0]);

      auto dx = func.Reverse(1, weight);
      expect_equal_eps(*d, dx[0], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[1], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[2], eps_deriv);
      d++;
    }
  }

  test_that("probit_integral gives the correct result"){
    /*
     integrand <- function(x, mu, sigma, k, s, use_log = FALSE){
    f1 <- dnorm(x, mean = mu, sd = sigma, log = use_log)
    f2 <- -pnorm(k + s * x, log.p = TRUE)
    if(use_log) f1 + log(f2) else f1 * f2
    }
    psi <- function(mu, sigma, k, s)
    integrate(
      integrand, mu = mu, sigma = sigma, k = k, s = s,
    lower = -Inf, upper = Inf, rel.tol = 1e-10, abs.tol = 0)$value
    # the approximation should work well for these values
    dput(mapply(
        psi,
    mu    =  c(0, -1, 1, 0, -1),
    sigma =  c(1, .5, 1, 2, .7),
    k     =  c(1, -1, 2, 1, .8),
    s     = -c(1, 1, 1, 1, 1)))
    dpsi <- function(mu, sigma, k){
    func <- function(x)
    psi(mu = x[1], sigma = x[2], k = x[3], s = -1)
    numDeriv::jacobian(func, c(mu, sigma, k), method.args = list(eps = 1e-10))
    }
    dput(mapply(dpsi,
    mu    =  c(0, -1, 1, 0, -1),
    sigma =  c(1, .5, 1, 2, .7),
    k     =  c(1, -1, 2, 1, .8)))
    */
    constexpr unsigned const n_nodes = 40L;

    std::vector<double> const
        mu     = { 0, -1, 1, 0, -1 },
        sigma  = { 1, .5, 1, 2, .7 },
        k      = { 1, -1, 2, 1, .8 },
        ex_res = { 0.35934760083045, 0.771877623464394, 0.35934760083045,
                   0.942178782481152, 0.078745703164216 },
       derivs = { 0.414713892660567, 0.377255144182391, -0.414715858841479,
                  0.825026811085921, 0.311729351150906, -0.825026811091284, 0.414715858841496,
                  0.377255144185024, -0.414715858840827, 0.686685373662083, 0.795828849179965,
                  -0.686690855788584, 0.133023446854648, 0.131235495730402, -0.133023446854061 };

    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = probit_integral(
        mu[i], sigma[i], k[i], n_nodes);
      expect_equal(ex_res[i], intval);
    }

    /* check Jacobian */
    using ADd = AD<double>;
    vector<ADd > a(3), b(1);
    a[0] = ADd(mu[0]);
    a[1] = ADd(sigma[0]);
    a[2] = ADd(k[0]);

    CppAD::Independent(a);
    b[0] = probit_integral(a[0], a[1], a[2], n_nodes);
    CppAD::ADFun<double> func(a, b);

    auto d = derivs.cbegin();
    vector<double> aa(3), weight(1);
    weight[0] = 1.;
    double const eps_deriv = std::pow(
      std::numeric_limits<double>::epsilon(), 1./ 4.);
    for(unsigned i = 0; i < mu.size(); ++i){
      aa[0] = mu[i];
      aa[1] = sigma[i];
      aa[2] = k[i];

      auto yy = func.Forward(0, aa);
      expect_equal(ex_res[i], yy[0]);

      auto dx = func.Reverse(1, weight);

      expect_equal_eps(*d, dx[0], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[1], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[2], eps_deriv);
      d++;
    }
  }
}
