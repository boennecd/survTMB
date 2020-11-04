#ifndef DMVNORM_LOG_H
#define DMVNORM_LOG_H

#include "utils.h"

namespace survTMB {
template<class Type>
Type mult_var_dens(vector<Type> const &theta, matrix<Type> const &rngs){
  if(theta.size() < 2L){ /* univariate case */
    Type const sd = exp(theta[0]),
      half(.5);
    Type out(0);
    Type const * const end = rngs.data() + rngs.size();
    for(Type const *x = rngs.data(); x != end; ++x){
      Type const scaled = *x / sd;
      out += - half * scaled * scaled;
    }

    out += - Type(rngs.size()) * log(Type(sqrt(2 * M_PI)) * sd);

    return out;
  }

  matrix<Type> const vcov = get_vcov_from_trian(theta);
  density::MVNORM_t<Type> norm_dist(vcov, true);
  Type out(0);
  for(unsigned g = 0; g < rngs.cols(); ++g)
    out -= norm_dist(rngs.col(g));

  return out;
}
} // namespace survTMB

#endif
