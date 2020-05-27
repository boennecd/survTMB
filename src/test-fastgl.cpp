#include "fastgl.h"
#include "testthat-wrap.h"
#include <limits>
#include "tmb_includes.h"

namespace {
/* int_-2^2 exp(x) dx */
template<class Type>
Type get_fastgl_testval(size_t const n){
  Type const ub( 2),
             lb(-2),
             T2(2.);
  auto const &xw = fastgl::GLPairsCached<Type>(n);

  Type out(0.);
  for(size_t i = 0; i < n; ++i){
    Type const node = (ub - lb) / T2 * xw[i].x + (ub + lb) / T2;
    out += exp(node) * xw[i].weight;
  }

  return out * (ub - lb) / T2;
}

} // namespace

context("Gauss-Legendre unit tests") {
  test_that("GLPairsCached yields correct value of a integral") {
    double const expect_val = exp(2) - exp(-2);

    for(unsigned n = 10L; n < 30L; ++n){
      double out = get_fastgl_testval<double>(n);
      expect_equal(out, expect_val);

      out = asDouble(get_fastgl_testval<CppAD::AD<double> >(n));
      expect_equal(out, expect_val);
    }
  }
}
