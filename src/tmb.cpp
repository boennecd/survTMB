#ifdef WITH_LIBTMB
#undef WITH_LIBTMB
#endif

/* we can either define all the atomic functions or only those that we need.
 * The latter results in a smaller .so object */
// #define TMB_PRECOMPILE

#include <TMB.hpp>

/* define just the atomic functions that we need */
#define PRECOMPILE_TMB_MACRO(ATOMIC_NAME)                      \
template                                                       \
  void ATOMIC_NAME<double>(const CppAD::vector<double>& tx,    \
                           CppAD::vector<double>& ty);         \
template                                                       \
  CppAD::vector<double> ATOMIC_NAME<double>                    \
  (const CppAD::vector<double>& tx)

namespace atomic {
PRECOMPILE_TMB_MACRO(pnorm1);
PRECOMPILE_TMB_MACRO(invpd);
PRECOMPILE_TMB_MACRO(matinv);
PRECOMPILE_TMB_MACRO(logdet);
PRECOMPILE_TMB_MACRO(matmul);

} // namespace atomic
