#include "testthat-wrap.h"
#include "utils.h"
#include <limits>
#include <vector>

using namespace survTMB;

context("Utils unit tests") {
  test_that("get_vcov_from_trian (MD) gives the correct result") {
    {
      /*
       Sigma <- matrix(c(4, 2, 1, 2, 2, .5, 1, .5, 1), 3L, 3L)
       ch <- t(chol(Sigma))
       diag(ch) <- log(diag(ch))
       dput(ch[lower.tri(ch, TRUE)])
       */
      vector<double> theta(6);
      theta << 0.693147180559945, 1, 0.5, 0, 0, -0.143841036225891;
      auto Sigma = get_vcov_from_trian(&theta[0L], 3L);

      std::vector<double> const ex { 4, 2, 1, 2, 2, .5, 1, .5, 1 };

      for(unsigned i = 0; i < ex.size(); ++i)
        expect_equal(ex[i], *(Sigma.data() + i));
    }

    {
      /*
       Sigma <- matrix(c(4, .33, .75, .33, 2, .5, .75, .5, 1), 3L, 3L)
       ch <- t(chol(Sigma))
       diag(ch) <- log(diag(ch))
       dput(ch <- ch[lower.tri(ch, TRUE)])

       fn <- function(x){
       o <- matrix(0, 3L, 3L)
       o[lower.tri(o, TRUE)] <- x
       diag(o) <- exp(diag(o))
       tcrossprod(o)
       }
       K <- matrixcalc::commutation.matrix(3L)
       dput(round(numDeriv::jacobian(fn, ch), 5))
       numDeriv::jacobian(fn, ch)
       */
      using ADd = AD<double>;
      constexpr size_t const n = 3L,
                         n_ele = (n * (n + 1L)) / 2L,
                            nn = n * n;
      vector<ADd> x(n_ele), y(nn);
      x[0] = ADd(0.693147180559945);
      x[1] = ADd(0.165);
      x[2] = ADd(0.375);
      x[3] = ADd(0.339720590501886);
      x[4] = ADd(0.31193151711395);
      x[5] = ADd(-0.13585598562127);

      CppAD::Independent(x);
      auto res = get_vcov_from_trian(x.data(), n);
      for(size_t i = 0; i < nn; ++i)
        *(y.data() + i) = *(res.data() + i);

      CppAD::ADFun<double> func(x, y);

      vector<double> xx(n_ele);
      for(size_t i = 0; i < n_ele; ++i)
        xx[i] = asDouble(x[i]);
      func.Forward(0, xx);

      double const eps_deriv = std::pow(
        std::numeric_limits<double>::epsilon(), 1./ 4.);
      std::vector<double> const expect =
              { 8, 0.33, 0.75, 0.33, 0, 0, 0.75, 0, 0, 0, 2, 0, 2,
                0.33, 0.375, 0, 0.375, 0, 0, 0, 2, 0, 0, 0.165, 2, 0.165, 0.75,
                0, 0, 0, 0, 3.94555, 0.43812, 0, 0.43812, 0, 0, 0, 0, 0, 0, 1.40456,
                0, 1.40456, 0.62386, 0, 0, 0, 0, 0, 0, 0, 0, 1.52415 };

      vector<double> w(nn);

#define CHECK_ELE_R(r)                                                \
      {                                                               \
        w.setZero();                                                  \
        w[r] = 1.;                                                    \
        auto dx = func.Reverse(1, w);                                 \
        expect_true(dx.size() == n_ele);                              \
        for(size_t j = 0; j < n_ele; ++j)                             \
          expect_equal_eps(expect[r + j * nn], dx[j], eps_deriv);     \
      }

     CHECK_ELE_R(0L)
     CHECK_ELE_R(1L)
     CHECK_ELE_R(2L)
     CHECK_ELE_R(3L)
     CHECK_ELE_R(4L)
     CHECK_ELE_R(5L)
     CHECK_ELE_R(6L)
     CHECK_ELE_R(7L)
     CHECK_ELE_R(8L)

#undef CHECK_ELE_R
    }

    {
      /*
       Z <- structure(c(0.94, 0.5, 0.19, -0.3, -0.42, 0.55, 0.89, 0.97, -0.7), .Dim = c(3L, 3L))
       dput(Z)
       Sigma <- matrix(c(4, .33, .75, .33, 2, .5, .75, .5, 1), 3L, 3L)
       ch <- t(chol(Sigma))
       diag(ch) <- log(diag(ch))
       dput(ch <- ch[lower.tri(ch, TRUE)])

       fn <- function(x){
       Z <- matrix(head(x, 9), 3, 3)
       o <- matrix(0, 3L, 3L)
       o[lower.tri(o, TRUE)] <- tail(x, -9)
       diag(o) <- exp(diag(o))
       O <- tcrossprod(o)
       sum(exp(Z * O))
       }
       K <- matrixcalc::commutation.matrix(3L)
       dput(round(numDeriv::jacobian(fn, c(Z, ch)), 8))
       numDeriv::jacobian(fn, c(Z, ch))

       dput(round(numDeriv::jacobian(fn, c(Z, ch)), 8))
       */
      using ADdd = AD<AD<double> >;
      using ADd =     AD<double>;
      constexpr size_t const n = 3L,
                         n_ele = (n * (n + 1L)) / 2L,
                            nn = n * n;
      std::vector<double> const vals_use =
        { 0.94, 0.5, 0.19, -0.3, -0.42, 0.55, 0.89, 0.97, -0.7,
          0.693147180559945, 0.165, 0.375, 0.339720590501886, 0.31193151711395,
          -0.13585598562127 };
      vector<ADdd> x(nn + n_ele), y(1L);

      /* Z */
      x.setZero();

      CppAD::Independent(x);
      matrix<ADdd> Z(n, n);
      {
        ADdd const * xp = &x[0];
        for(size_t c = 0; c < n; ++c)
          for(size_t r = 0; r < n; ++r, ++xp)
            Z(r, c) = *xp;
      }

      auto S = get_vcov_from_trian(x.data() + nn, n);
      y[0L] = ADdd(0);
      for(size_t c = 0; c < n; ++c)
        for(size_t r = 0; r < n; ++r)
          y[0L] += exp(Z(r, c) * S(r, c));

      CppAD::ADFun<ADd> func_ADdd;
      func_ADdd.Dependent(x, y);
      func_ADdd.optimize();

      vector<ADd> xx(x.size());
      for(unsigned i = 0; i < x.size(); ++i)
        xx[i] = CppAD::Value(x[i]);

      CppAD::Independent(xx);
      vector<ADd> yy = func_ADdd.Jacobian(xx);

      CppAD::ADFun<double> func_ADd;
      func_ADd.Dependent(xx, yy);
      func_ADd.optimize();

      // check the gradient
      double const eps_deriv = std::pow(
        std::numeric_limits<double>::epsilon(), 1./ 4.);
      {
        std::vector<double> const expect =
          { 171.79370392, 0.38919973, 0.86486481, 0.29889509,
            0.86342105, 0.65826534, 1.46201837, 0.8120875, 0.4965853, 324.54261536,
            1.43844052, 4.02677219, 0.29208579, 3.01297205, -0.52980846 };

        vector<ADd> xx_use(vals_use.size());
        for(size_t i = 0; i < vals_use.size(); ++i)
          xx_use[i] = ADd(vals_use[i]);

        auto dx = func_ADdd.Jacobian(xx_use);
        expect_true(static_cast<size_t>(dx.size()) == expect.size());
        for(size_t i = 0; i < nn + n_ele; ++i)
          expect_equal_eps(expect[i], CppAD::Value(dx[i]), eps_deriv);
      }

      // check the Hessian
      {
        std::vector<double> const expect =
          { 687.17481566, 0, 0, 0, 0, 0, 0, 0, 0, 1635.47603406,
            0, 0, 0, 0, 0, 0, 0.12843591, 0, 0, 0, 0, 0, 0, 0, 0.45341768,
            2.74798597, 0, 0, 0, 0, 0, 0, 0.64864861, 0, 0, 0, 0, 0, 0, 0.98810805,
            0, 2.63495479, 0, 0, 0, 0, 0, 0, 0.09863538, 0, 0, 0, 0, 0, 0.26930448,
            1.63214836, 0, 0, 0, 0, 0, 0, 0, 0, 1.72684209, 0, 0, 0, 0, 0,
            0.02279432, 0, 0.27253367, 0, 0, 0, 0, 0, 0, 0, 0.32913267, 0,
            0, 0, 0, 0.62946623, 0.27696514, 0.73542638, 2.35765332, 0, 0,
            0, 0, 0, 0, 0, 1.09651378, 0, 0, 2.43791564, 0, 6.50110836, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0.40604375, 0, 0, 0.90446246, 0.39796348,
            1.05671364, 3.38764626, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4965853,
            0, 0, 0.11173169, 0, 0.09294036, 0.22706077, 1635.47603406, 0.45341768,
            0.98810805, 0.26930448, 0, 0, 2.43791564, 0, 0, 3077.19839748,
            0.88434848, 6.28662785, 0, 0, 0, 0, 2.74798597, 0, 1.63214836,
            0.02279432, 0.62946623, 0, 0.90446246, 0, 0.88434848, 1.42202197,
            2.41873991, 0.41566211, 1.01466998, 0, 0, 0, 2.63495479, 0, 0,
            0.27696514, 6.50110836, 0.39796348, 0.11173169, 6.28662785, 2.41873991,
            5.83695979, 0.13926332, 0.56030674, 0.27814944, 0, 0, 0, 0, 0.27253367,
            0.73542638, 0, 1.05671364, 0, 0, 0.41566211, 0.13926332, 1.1319842,
            4.41530566, 0, 0, 0, 0, 0, 0, 2.35765332, 0, 3.38764626, 0.09294036,
            0, 1.01466998, 0.56030674, 4.41530566, 3.19991095, 0.23136954,
            0, 0, 0, 0, 0, 0, 0, 0, 0.22706077, 0, 0, 0.27814944, 0, 0.23136954,
            -0.49436257 };

        vector<double> xxx(xx.size());
        for(unsigned i = 0; i < xx.size(); ++i)
          xxx[i] = vals_use[i];

        vector<double> hes = func_ADd.Jacobian(xxx);
        expect_true(static_cast<size_t>(hes.size()) == expect.size());
        for(size_t c = 0; c < static_cast<size_t>(x.size()); ++c)
          for(size_t r = 0; r < static_cast<size_t>(x.size()); ++r)
            expect_equal_eps(expect[r + c * x.size()],
                             hes[r + c * x.size()], eps_deriv);
      }
    }
  }

  test_that("get_vcov_from_trian (1D) gives the correct result") {
    double theta = -1;
    auto Sigma = get_vcov_from_trian(&theta, 1L);

    double ex = std::exp(2 * theta);
    expect_equal(ex, *Sigma.data());
  }
}
