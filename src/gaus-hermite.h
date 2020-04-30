#ifndef GAUS_HERMITE_H
#define GAUS_HERMITE_H
#include <vector>
#include "tmb_includes.h"

namespace GaussHermite {

template<class Type = double>
struct HermiteData {
  std::vector<Type> x, w;

  HermiteData(unsigned const n): x(n), w(n) { }

  template<class O_Type>
  HermiteData(HermiteData<O_Type> const &o){
    std::size_t const n = o.x.size();
    x.resize(n);
    w.resize(n);

    for(unsigned i = 0; i < n; ++i){
      x[i] = asDouble(o.x[i]);
      w[i] = asDouble(o.w[i]);
    }
  }
};

HermiteData<double> GaussHermiteData(unsigned const);

template<class Type>
HermiteData<Type> const& GaussHermiteDataCached(unsigned const);

constexpr std::size_t GaussHermiteDataCachedMaxArg(){
  return 100L;
}

} // namespace GaussHermite

#endif
