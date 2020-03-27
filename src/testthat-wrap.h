#ifndef TESTTHAT_WRAP_H
#define TESTTHAT_WRAP_H
#include <testthat.h>
#include <limits>

template<class T1, class T2>
void expect_equal(T1 const target, T2 const current){
  double const eps = std::sqrt(std::numeric_limits<double>::epsilon());
  if(std::abs(target) < eps)
    expect_true(std::abs( target - current          ) < eps);
  else
    expect_true(std::abs((target - current) / target) < eps);
}

#endif
