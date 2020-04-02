#ifndef LAPLACE_H
#define LAPLACE_H

#include "tmb_includes.h"
#include "common.h"

namespace survTMB {
/* Computes the log-likelihood for given random effects.
 *
 * Args:
 *   COMMON_ARGS: see the other header file.
 *   u: [random effect dim] x [# groups] matrix with random effects.
 */
template<class Type>
void laplace(COMMON_ARGS(Type, parallel_accumulator), matrix<Type> const &u);

} // namespace survTMB

#endif
