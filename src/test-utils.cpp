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
  }

  test_that("get_vcov_from_trian (1D) gives the correct result") {
    double theta = -1;
    auto Sigma = get_vcov_from_trian(&theta, 1L);

    double ex = std::exp(2 * theta);
    expect_equal(ex, *Sigma.data());
  }
}
