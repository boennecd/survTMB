#ifndef GVA_H
#define GVA_H

#include "gaus-hermite.h"
#include "common.h"

namespace survTMB {

/* Computes the lower bound.
 *
 * Args:
 *   COMMON_ARGS: see the other header file.
 *   theta_VA: vector with VA parameters for each group in terms. For each
 *             group, the first elements are the mean and the subsequent
 *             elemenet is an upper diagonal matrix L such that the
 *             covariance matrix is L^\top L.
 */
template<class Type, template <class> class Accumlator>
void GVA(COMMON_ARGS(Type, Accumlator), vector<Type> const &theta_VA,
         unsigned const n_nodes);

} // namespace survTMB

#endif
