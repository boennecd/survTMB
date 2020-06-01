#ifndef GAMMA_TO_NU_H
#define GAMMA_TO_NU_H

#include "tmb_includes.h"

/* quantity needed when mapping from CP to DP w/ SNVA. R functions are
func <- function(g){
g_abs <- abs(g)
cv <- 2 * g_abs / (4 - pi)
out <-  cv^(1/3) / sqrt(1 + cv^(2/3))
ifelse(g < 0, -out, out)
}
dfunc <- function(g){
cv <- 2 * g / (4 - pi)
cv_2_3 <- (cv * cv)^(1/3)
1 / (3 * (cv_2_3 + 1)^(3/2) * cv_2_3) * 2 / (4 - pi)
}
plot(func, xlim = c(-.9, .9))
gs <- seq(-.9, .9, length.out = 100)
all.equal(sapply(gs, numDeriv::grad, func = func),
          dfunc(gs))
*/
inline double gamma_to_nu_func(double const g){
  constexpr double mult = 2. / (4. - M_PI);
  double const g_sign = g < 0 ? -1. : 1.,
               g_abs  =  g * g_sign,
               cv     = mult * g_abs,
               out    = pow(cv, 1. / 3.) / sqrt(1. + pow(cv, 2. / 3.));
  return g_sign * out;
}

template<class Type>
Type dgamma_to_nu_func(Type const g){
  Type const mult(2. / (4. - M_PI)),
             cv = mult * g,
         cv_2_3 = pow(cv * cv, 1. / 3.),
          denom = 3. * pow(cv_2_3 + 1., 3. / 2.) * cv_2_3,
            out = mult / denom,
            eps = Type(1. / std::numeric_limits<double>::epsilon());

  return CppAD::CondExpLe(out, eps, out, eps);
}

namespace atomic {
TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  gamma_to_nu1
  ,
  // OUTPUT_DIM
  1
  ,
  // ATOMIC_DOUBLE
  ty[0] =  gamma_to_nu_func(tx[0]);
  ,
  // ATOMIC_REVERSE
  px[0] = dgamma_to_nu_func(tx[0]) * py[0];
  )

} // namespace atomic

template<class Type>
Type gamma_to_nu(Type g){
  CppAD::vector<Type> tx(1);
  tx[0] = g;
  return atomic::gamma_to_nu1(tx)[0];
}

VECTORIZE1_t(gamma_to_nu)

#endif
