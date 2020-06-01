#include "testthat-wrap.h"
#include "snva-utils.h"
#include <vector>

using namespace GaussHermite;
using namespace GaussHermite::SNVA;

namespace {
template<class C1, class C2>
void do_test(C1 const &ex, C2 const &re){
  expect_true(ex.size() == (unsigned)re.size());
  for(unsigned i = 0; i < ex.size(); ++i)
    expect_equal(ex[i], *(re.data() + i));
}
} // namespace

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
     dput(psi(c(0:2, 100)))
     library(numDeriv)
     dput(sapply(c(0:2, 100), function(x)
     grad(psi, x)))
     */
    constexpr size_t const n_tests(4L);
    constexpr double const sigmas[n_tests] = { 0., 1., 2., 100. };
    constexpr double const expect[n_tests] = {
      -0.693147180559945, -0.5, -0.319885055354923, -0.00720608477031823 };
    constexpr double const dx    [n_tests] = {
      33690.3759242532, 0.24026897903233, 0.126022609306171,
      7.20531098254606e-05 };
    constexpr unsigned const n_nodes = 20L;

    using ADd = AD<double>;
    vector<ADd > a(1), b(1);
    a[0] = ADd(1.5);
    CppAD::Independent(a);
    a[0] = a[0] * a[0];
    b[0] = entropy_term(a[0], n_nodes);
    CppAD::ADFun<double> func(a, b);

    vector<double> x(1), w(1);
    w[0] = 1;

    for(unsigned i = 0L; i < n_tests; ++i){
      double res = entropy_term(sigmas[i] * sigmas[i], n_nodes);
      expect_equal(expect[i], res);

      x[0] = sigmas[i];
      auto const y = func.Forward(0, x);
      expect_equal(expect[i], y[0]);

      if(i > 0){
        auto dy = func.Reverse(1, w);
        expect_equal(dx[i], dy[0]);
      }
    }

    {
      double res = entropy_term(1e-8, n_nodes);
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
     dpsi <- function(mu, sigma, rho, k){
     func <- function(x)
     psi(mu = x[1], sigma = x[2], rho = x[3], k = x[4], s = 1)
     numDeriv::jacobian(func, c(mu, sigma, rho, k),
     method.args = list(eps = 1e-10))
     }
     dput(mapply(dpsi,
     mu    = c(0, -1, 1, 0, -1),
     sigma = c(1, .5, 1, 2, .7),
     rho   = c(0, -1, 1, 1, -2),
     k     = c(1, .3, 2, 1, .8)))
     */
    constexpr unsigned const n_nodes = 40L;
    std::vector<double> const
        mu     = { 0, -1, 1, 0, -1 },
        sigma  = { 1, .5, 1, 2, .7 },
        rho    = { 0, -1, 1, 1, -2 },
        k      = { 1, .3, 2, 1, .8 },
        ex_res = { 0.80605918334744, 0.0969646948660297, 2.38729257810521, 1.78136403451249,
                   0.189741259960398 },
        derivs = { 0.500006490454055, 0.206620964134777, 0.398937704544547,
                   0.499999999989764, 0.0915328990478164, -0.0276112082285875, 0.0152040703474796,
                   0.305109663495911, 0.882551315240342, 0.793854072803876, 0.232363697133678,
                   0.44127565761652, 0.73931997580623, 0.731042542671959, 0.142729929312629,
                   0.739329245111041, 0.169614547405035, -0.10904339875669, 0.0180383959195121,
                   0.212018184253227 };

    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = mlogit_integral(
        mu[i], sigma[i], rho[i], std::log(k[i]), n_nodes);
      expect_equal(ex_res[i], intval);
    }

    /* check Jacobian */
    using ADd = AD<double>;
    vector<ADd > a(4), b(1);
    a[0] = ADd(mu[0]);
    a[1] = ADd(sigma[0]);
    a[2] = ADd(rho[0]);
    a[3] = ADd(k[0]);

    CppAD::Independent(a);
    b[0] = mlogit_integral(a[0], a[1], a[2], log(a[3]), n_nodes);
    CppAD::ADFun<double> func(a, b);

    auto d = derivs.cbegin();
    vector<double> aa(4), weight(1);
    weight[0] = 1.;
    double const eps_deriv = std::pow(
      std::numeric_limits<double>::epsilon(), 1./ 4.);
    for(unsigned i = 0; i < mu.size(); ++i){
      aa[0] = mu[i];
      aa[1] = sigma[i];
      aa[2] = rho[i];
      aa[3] = k[i];

      auto yy = func.Forward(0, aa);
      expect_equal(ex_res[i], yy[0]);

      auto dx = func.Reverse(1, weight);

      expect_equal_eps(*d, dx[0], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[1], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[2], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[3], eps_deriv);
      d++;
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
     dpsi <- function(mu, sigma, rho, k){
     func <- function(x)
     psi(mu = x[1], sigma = x[2], rho = x[3], k = x[4], s = -1)
     numDeriv::jacobian(func, c(mu, sigma, rho, k),
     method.args = list(eps = 1e-10))
     }
     dput(mapply(dpsi,
     mu    =  c(0, -1,  1,  0,  -1),
     sigma =  c(1, .5,  1,  2,  .7),
     rho   =  c(0, -1,  1,  1,  -2),
     k     =  c(1, .3, -2,  1, -.8)))
     */
    unsigned const n_nodes = 40L;
    std::vector<double> const
        mu     = { 0, -1,  1,  0,  -1 },
        sigma  = { 1, .5,  1,  2,  .7 },
        rho    = { 0, -1,  1,  1,  -2 },
        k      = { 1, .3, -2,  1, -.8 },
        ex_res = { 0.35934760083045, 0.0981790817582183, 8.93028198953091, 1.79323321188706,
                   0.360643279313754 },
        derivs = { 0.414713892660567, 0.377255144182391, 0.330895837479333,
                   -0.414715858841479, 0.171778814619194, -0.0170681853114501, 0.0309884723079667,
                   -0.171778814619795, 3.82006109597964, 3.72538370057417, 0.928490173628573,
                   -3.82006109586034, 1.24521443636393, 1.62553040226332, 0.111857457341852,
                   -1.24520007126255, 0.471920179688172, -0.278392484375124, 0.0533432871738786,
                   -0.471920179691003 };
    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = probit_integral(
        mu[i], sigma[i], rho[i], k[i], n_nodes);

      expect_equal(ex_res[i], intval);
    }

    /* check Jacobian */
    using ADd = AD<double>;
    vector<ADd > a(4), b(1);
    a[0] = ADd(mu[0]);
    a[1] = ADd(sigma[0]);
    a[2] = ADd(rho[0]);
    a[3] = ADd(k[0]);

    CppAD::Independent(a);
    b[0] = probit_integral(a[0], a[1], a[2], a[3], n_nodes);
    CppAD::ADFun<double> func(a, b);

    auto d = derivs.cbegin();
    vector<double> aa(4), weight(1);
    weight[0] = 1.;
    double const eps_deriv = std::pow(
      std::numeric_limits<double>::epsilon(), 1./ 4.);
    for(unsigned i = 0; i < mu.size(); ++i){
      aa[0] = mu[i];
      aa[1] = sigma[i];
      aa[2] = rho[i];
      aa[3] = k[i];

      auto yy = func.Forward(0, aa);
      expect_equal(ex_res[i], yy[0]);

      auto dx = func.Reverse(1, weight);

      expect_equal_eps(*d, dx[0], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[1], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[2], eps_deriv);
      d++;
      expect_equal_eps(*d, dx[3], eps_deriv);
      d++;
    }
  }

  test_that("SNVA_MD_theta_DP_to_DP gives the correct result") {
    /*
     dput(mu <- c(4, 3, 8))
     dput(Sig <- matrix(c(1, 1, 2, 1, 4, 3, 2, 3, 6), 3))
     rho <- c(5, 7, 9)
     dput(alpha <- rho * sqrt(diag(Sig)))
     ch <- t(chol(Sig))
     log_sigs <- log(diag(ch))
     L <- diag(diag(ch)^(-1)) %*% ch
     dput(c(mu, log_sigs, L[lower.tri(L)], rho))
     */
    using Type = double;
    vector<Type> theta(12);
    theta << 4, 3, 8, 0, 0.549306144334055, 0.255412811882995,
             0.577350269189626,
             1.54919333848297, 0.447213595499958,
             5, 14, 22.0454076850486;

    auto input = SNVA_MD_theta_DP_to_DP(&theta[0], theta.size(), 3L);
    std::vector<double> const mu = { 4, 3, 8 },
                             Sig = { 1, 1, 2, 1, 4, 3, 2, 3, 6 },
                             rho = { 5, 7, 9 };

    expect_true(input.va_mus.size() == 1L);

    do_test(mu , input.va_mus[0]);
    do_test(Sig, input.va_lambdas[0]);
    do_test(rho, input.va_rhos[0]);
  }
  {
    /*
     Sigu <- c(4, 2, 1, 1, 2, 3)
     Sig <- matrix(0., 3L, 3L)
     Sig[upper.tri(Sig, TRUE)] <- Sigu
     Sig <- crossprod(Sig)
     gamma_trans <- c(0.01, -0.01, 0)
     skew_boundary <- 0.99527
     gam <- 2 * skew_boundary / (1 + exp(-gamma_trans)) - skew_boundary

     dput(cov_to_theta(Sig),  control = c("keepNA", "keepInteger"))
     dput(out <- survTMB:::cp_to_dp(
     mu = c(-1, 0, 1),
     Sigma = Sig,
     gamma = gam))
     dput(out$alpha / sqrt(diag(out$Psi)))
     */

    using Type = double;
    vector<Type> theta(12);
    theta << -1., 0., 1.,
             1.38629436111989, 0, 1.09861228866811, 2, 0.333333333333333,
             0.666666666666667,
             0.01, -0.01, 0;

    auto input = SNVA_MD_theta_CP_trans_to_DP(&theta[0], theta.size(), 3L);
    std::vector<double> const
      mu  = { -1.90533215899427, 0.506096062431957, 1 },
      Sig = { 16.8196263181092,
              7.54181495913998, 4, 7.54181495913998, 5.25613322440913, 4, 4,
              4, 14 },
      rho = { 1.09012077024187, -2.14930211613867,
              0.302623241684799 };

    expect_true(input.va_mus.size() == 1L);

    do_test(mu , input.va_mus[0]);
    do_test(Sig, input.va_lambdas[0]);
    do_test(rho, input.va_rhos[0]);
  }
}
