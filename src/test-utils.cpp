#include <testthat.h>
#include "utils.h"
#include <limits>
#include <vector>

using namespace survTMB;

context("Utils unit tests") {
  test_that("get_vcov_from_trian (MD) gives the correct result") {
    /*
     Sigma <- matrix(c(4, 2, 1, 2, 2, .5, 1, .5, 1), 3L, 3L)
     ch <- t(chol(Sigma))
     sigs <- diag(ch)
     L <- diag(sigs^(-1)) %*% ch
     dput(c(log(sigs), L[lower.tri(L)]))
     */
    vector<double> theta(6);
    theta << 0.693147180559945, 0, -0.143841036225891, 1,
             0.577350269189626, 0;
    auto Sigma = get_vcov_from_trian(&theta[0L], 3L);

    std::vector<double> const ex { 4, 2, 1, 2, 2, .5, 1, .5, 1 };
    constexpr double eps = std::sqrt(
      std::numeric_limits<double>::epsilon());

    for(unsigned i = 0; i < ex.size(); ++i)
      expect_true(std::abs((ex[i] - *(Sigma.data() + i)) / ex[i]) < eps);
  }
  test_that("get_vcov_from_trian (1D) gives the correct result") {
    double theta = -1;
    auto Sigma = get_vcov_from_trian(&theta, 1L);

    double ex = std::exp(2 * theta);
    constexpr double eps = std::sqrt(
      std::numeric_limits<double>::epsilon());
    expect_true(std::abs((ex - *Sigma.data()) / ex) < eps);
  }
}
