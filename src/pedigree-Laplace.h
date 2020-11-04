#ifndef PEDIGREE_LAPLACE_H
#define PEDIGREE_LAPLACE_H

#include "tmb_includes.h"

namespace survTMB {

template<class Type>
void pedigree_laplace
  (parallel_accumulator<Type> &out, SEXP dat,
   vector<Type> const &omega, vector<Type> const &beta,
   vector<Type> const &log_sds, vector<Type> const &rng_modes);
}

#endif
