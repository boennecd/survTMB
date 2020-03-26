#ifndef TMB_INCLUDES_H
#define TMB_INCLUDES_H

#define WITH_LIBTMB

#ifdef TMB_PRECOMPILE
#undef TMB_PRECOMPILE
#endif

#ifdef CSKIP
#undef CSKIP
#endif

#include <TMB.hpp>

namespace CppAD {
inline double Value(double const x){
  return x;
}
}

#endif
