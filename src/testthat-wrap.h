#ifndef TESTTHAT_WRAP_H
#define TESTTHAT_WRAP_H
#include <testthat.h>
#include <limits>


#define expect_equal(TARGET, CURRENT)                          \
  {                                                            \
    double const eps =                                         \
      std::sqrt(std::numeric_limits<double>::epsilon());       \
                                                               \
    if(std::abs((TARGET)) > eps)                               \
      expect_true(std::abs((TARGET) - (CURRENT)) /             \
                    std::abs((TARGET)) < eps);                 \
    else                                                       \
      expect_true(std::abs((TARGET) - (CURRENT) < eps));       \
  }

#define expect_equal_eps(TARGET, CURRENT, EPS)                 \
{                                                              \
  if(std::abs((TARGET)) > (EPS))                               \
    expect_true(std::abs((TARGET) - (CURRENT)) /               \
      std::abs((TARGET)) < (EPS));                             \
  else                                                         \
    expect_true(std::abs((TARGET) - (CURRENT) < (EPS)));       \
}

#endif
