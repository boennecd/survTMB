#ifndef PNORM_LOG_H
#define PNORM_LOG_H

#include "tmb_includes.h"

namespace atomic {
/* Computes log CDF of standard normal distribution */
TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  pnorm_log1
  ,
  // OUTPUT_DIM
  1
  ,
  // ATOMIC_DOUBLE
  ty[0] = Rmath::Rf_pnorm5(tx[0], 0, 1, 1, 1);
  ,
  // ATOMIC_REVERSE
  px[0] = py[0];
  if(CppAD::Value(tx[0]) > -10){
    Type const cdf = exp(ty[0]);
    px[0] *= dnorm1(tx[0]) / cdf;
  } else {
    Type const log_pdf = dnorm(tx[0], Type(0.), Type(1.), 1L);
    px[0] *= exp(log_pdf - ty[0]);
  })
} // namespace atomic

/* Computes the log CDF of normal distribution.
 *
 * Args:
 *   Similar to `stats::pnorm`.
 */
template<class Type>
Type pnorm_log(Type q, Type  mean = 0., Type sd = 1.){
  CppAD::vector<Type> tx(1);
  tx[0] = (q - mean) / sd;
  return atomic::pnorm_log1(tx)[0];
}
VECTORIZE3_ttt(pnorm_log)
VECTORIZE1_t  (pnorm_log)

#endif
