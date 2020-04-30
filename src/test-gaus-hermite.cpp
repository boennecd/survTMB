#include "gaus-hermite.h"
#include "testthat-wrap.h"
#include <limits>

context("gaus-hermite unit tests") {
  test_that("GaussHermiteData gives the correct result") {
    /*
     xw <- fastGHQuad::gaussHermiteData(5)
     paste("{",
     paste0(sprintf("%20.16f", xw$x, xw$x), collapse = ", "),
     "}")
     paste("{",
     paste0(sprintf("%20.16f", xw$w, xw$x), collapse = ", "),
     "}")
     */
    constexpr unsigned n = 5L;
    auto xw = GaussHermite::GaussHermiteData(n);
    GaussHermite::HermiteData<CppAD::AD<double> > xw_AD(xw);
    std::vector<double> const x =
      {  -2.0201828704560856,  -0.9585724646138187,   0.0000000000000000,   0.9585724646138187,   2.0201828704560856 };
    std::vector<double> const w =
      {   0.0199532420590459,   0.3936193231522411,   0.9453087204829413,   0.3936193231522414,   0.0199532420590459 };

    expect_true(x.size() == xw.x.size());
    expect_true(w.size() == xw.w.size());

    for(unsigned i = 0; i < n; ++i){
      expect_equal(x[i], xw.x[i]);
      expect_equal(w[i], xw.w[i]);
      expect_equal(x[i], asDouble(xw_AD.x[i]));
      expect_equal(w[i], asDouble(xw_AD.w[i]));
    }
  }

  test_that("GaussHermiteData yields correct value of a integral") {
    /* int x phi(x, 2, 1) dx =  int (2 + sqrt(2) x) phi(x, 0, 2) dx = 2*/
    for(unsigned n = 10L; n < 30L; ++n){
      auto xw = GaussHermite::GaussHermiteData(n);
      double out(0.);
      for(unsigned i = 0; i < n; ++i)
        out += xw.w[i] * (2 + sqrt(2) * xw.x[i]);
      out /= sqrt(M_PI);

      expect_equal(out, 2.);
    }
  }
}
