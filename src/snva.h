#ifndef SNVA_H
#define SNVA_H

#include "gaus-hermite.h"
#include "common.h"

namespace survTMB {

/* Computes the lower bound.
 *
 * Args:
 *   COMMON_ARGS: see the other header file.
 *   TODO: document other parameters.
 */
template<class Type, template <class> class Accumlator>
void SNVA(COMMON_ARGS(Type, Accumlator), vector<Type> const &theta_VA,
          unsigned const n_nodes, std::string const &param_type);

} // namespace survTMB

#endif
