#include "testthat-wrap.h"
#include "gsm.h"
#include <vector>

namespace {
template<class Fam>
struct test_gsm_fam_util {
  std::vector<double> g_log, gp, gp_g, gpp_gp, gpp;

  test_gsm_fam_util(std::vector<double> const &eta){
    size_t const n = eta.size();
    g_log .reserve(n);
    gp    .reserve(n);
    gp_g  .reserve(n);
    gpp_gp.reserve(n);
    gpp   .reserve(n);

    for(auto &e : eta){
      Fam f(e);
      g_log .emplace_back(f.g_log());
      gp    .emplace_back(f.gp());
      gp_g  .emplace_back(f.gp_g());
      gpp_gp.emplace_back(f.gpp_gp());
      gpp   .emplace_back(f.gpp());
    }
  }
};
} // namespace

context("testing gsm") {
  test_that("ph link is correct") {
    /*
     dput(eta <- as.numeric((-3):3))
     g <- function(x)     exp(-exp(x))
     g_log <- function(x)     -exp(x)
     library(numDeriv)
     dput(g_log(eta))
     dput(gp <- sapply(eta, function(x) grad(g, x, method.args=list(eps = 1e-8))))
     dput(gpp <- sapply(eta, function(x) hessian(g, x, method.args=list(eps = 1e-8))))
     dput(gp / g(eta))
     dput(gpp / gp)
     */
    std::vector<double> const
      eta    = { -3, -2, -1, 0, 1, 2, 3 },
      g_log  = { -0.0497870683678639, -0.135335283236613, -0.367879441171442,
                 -1, -2.71828182845905, -7.38905609893065, -20.0855369231877 },
      gp     = { -0.0473690096769045, -0.118204951592349, -0.254646380044638,
                 -0.367879441170156, -0.179374078732767, -0.00456628142010774,
                 -3.80054250403406e-08 },
      gpp    = { -0.0450106455545372, -0.102207650989543, -0.160967212057699,
                 0, 0.308215219985128, 0.029174228156663, 7.25350521828551e-07 },
      gp_g   = { -0.0497870683668093, -0.135335283235704, -0.367879441172967,
                 -0.999999999996503, -2.7182818284401, -7.389056098898, -20.0855369231333 },
      gpp_gp = { 0.950212931651871, 0.86466471677112, 0.632120558829393, -1.50795421548183e-06,
                 -1.71828182847038, -6.38905609895034, -19.0854469081357 };

    test_gsm_fam_util<gsm_objs::gsm_ph> fam(eta);
    for(size_t i = 0; i < eta.size(); ++i){
      expect_equal_eps(g_log [i], fam.g_log [i], 1e-5);
      expect_equal_eps(gp    [i], fam.gp    [i], 1e-5);
      expect_equal_eps(gpp   [i], fam.gpp   [i], 1e-5);
      expect_equal_eps(gp_g  [i], fam.gp_g  [i], 1e-5);
      expect_equal_eps(gpp_gp[i], fam.gpp_gp[i], 1e-5);
    }
  }

  test_that("po link is correct") {
    /*
     dput(eta <- c(-40, (-3):3, 40))
     g <- function(x)     1 / (1 + exp(x))
     g_log <- function(x) -log(1 + exp(x))
     library(numDeriv)
     dput(g_log(eta))
     dput(gp <- sapply(eta, function(x) grad(g, x, method.args=list(eps = 1e-8))))
     dput(gpp <- sapply(eta, function(x) hessian(g, x, method.args=list(eps = 1e-8))))
     dput(gp / g(eta))
     dput(gpp / gp)
     */
    std::vector<double> const
    eta    = { -40, -3, -2, -1, 0, 1, 2, 3, 40 },
    g_log  = { 0, -0.048587351573742, -0.126928011042973, -0.313261687518223,
               -0.693147180559945, -1.31326168751822, -2.12692801104297, -3.04858735157374,
               -40 },
    gp     = { 0, -0.0451766597314526, -0.104993585403967, -0.196611933248946,
               -0.250000000003241, -0.196611933242019, -0.104993585403084, -0.0451766597307574,
               0 },
    gpp    = { -4.89516324790633e-21, -0.0408915746611961, -0.0799625010557186,
               -0.0908577476718339, 0, 0.0908577476729026, 0.0799625010561773,
               0.040891574660934, 4.89516324790633e-21 },
    gp_g   = { 0, -0.0474258731781341, -0.11920292202264, -0.268941421380205,
               -0.500000000006481, -0.731058578632004, -0.880797077974339, -0.952574126819171,
               -1.00000000000779 },
    gpp_gp = { 1, 0.90514825363963, 0.761594155948289, 0.462117157236798,
               -3.36849889541162e-07, -0.462117157258513, -0.761594155959059,
               -0.905148253647757, -1 };

    test_gsm_fam_util<gsm_objs::gsm_logit> fam(eta);
    for(size_t i = 0; i < eta.size(); ++i){
      expect_equal_eps(g_log [i], fam.g_log [i], 1e-5);
      expect_equal_eps(gp    [i], fam.gp    [i], 1e-5);
      expect_equal_eps(gpp   [i], fam.gpp   [i], 1e-5);
      expect_equal_eps(gp_g  [i], fam.gp_g  [i], 1e-5);
      expect_equal_eps(gpp_gp[i], fam.gpp_gp[i], 1e-5);
    }
  }

  test_that("probit link is correct") {
    /*
     dput(eta <- c(-40, (-3):3, 40))
     g <- function(x)     pnorm(-x)
     g_log <- function(x) pnorm(-x, log.p = TRUE)
     library(numDeriv)
     dput(g_log(eta))
     dput(gp <- sapply(eta, function(x) grad(g, x, method.args=list(eps = 1e-8))))
     dput(gpp <- sapply(eta, function(x) hessian(g, x, method.args=list(eps = 1e-8))))
     dput(gp / g(eta))
     dput(gpp / gp)
     */
    std::vector<double> const
    eta    = { -40, -3, -2, -1, 0, 1, 2, 3, 40 },
    g_log  = { 0, -0.00135080996474819, -0.0230129093289635, -0.17275377902345,
               -0.693147180559945, -1.84102164500926, -3.78318433368203, -6.60772622151035,
               -804.608442013754 },
    gp     = { 0, -0.00443184841168233, -0.0539909665111145, -0.241970724519369,
               -0.398942306029329, -0.24197072451802, -0.0539909665128562, -0.00443184841192828,
               0 },
    gpp    = { 0, -0.0132955452357521, -0.107981933026381, -0.241970724519003,
               0, 0.241970724519496, 0.10798193302636, 0.0132955452357508, -9.22095252600812e-289 },
    gp_g   = { 0, -0.00443783904186964, -0.0552478626768681, -0.287599970939447,
               -0.797884612058659, -1.5251352761539, -2.37321553280826, -3.28309865492323,
               -40.0249688472063 },
    gpp_gp = { 40., 3.00000000015909, 2.0000000000769, 0.999999999998486,
               0, -1.0000000000061, -2.00000000001199, -2.99999999999231,
               -40 };

    test_gsm_fam_util<gsm_objs::gsm_probit> fam(eta);
    for(size_t i = 0; i < eta.size(); ++i){
       expect_equal_eps(g_log [i], fam.g_log [i], 1e-5);
       expect_equal_eps(gp    [i], fam.gp    [i], 1e-5);
       expect_equal_eps(gpp   [i], fam.gpp   [i], 1e-5);
       expect_equal_eps(gp_g  [i], fam.gp_g  [i], 1e-5);
       expect_equal_eps(gpp_gp[i], fam.gpp_gp[i], 1e-5);
    }
  }
}
