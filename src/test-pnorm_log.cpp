#include "testthat-wrap.h"
#include "pnorm-log.h"
#include <vector>
#include <limits>

context("pnorm_log unit tests") {
  test_that("pnorm_log gives the correct result") {
    using std::abs;
    using std::vector;
    using CppAD::AD;

    /*
     x <- 2; mu <- 1; sd <- 2
     dput(pnorm(x, mu, sd, log.p = TRUE))
     numericDeriv(quote(pnorm(x, mu, sd, log.p = TRUE)), theta = "x")
     dput(exp(dnorm(x, mu, sd, 1) - pnorm(x, mu, sd, log.p = TRUE)))
     */
    vector<AD<double> > x(1);
    x[0] = 3;
    Independent(x);
    AD<double> const mu = 1, sd = 2;
    vector<AD<double> > y(1);
    y[0] = pnorm_log(x[0], mu, sd);
    CppAD::ADFun<double> func(x, y);

    vector<double> xx(1);
    xx[0] = 2;
    vector<double> yy = func.Forward(0, xx);
    constexpr double const true_yy(-0.368946415288657),
                           true_dy(0.254580216918517);
    expect_equal(true_yy, yy[0]);

    vector<double> w(1);
    w[0] = 1;
    vector<double> dy = func.Reverse(1, w);
    expect_equal(true_dy, dy[0]);
  }
}
