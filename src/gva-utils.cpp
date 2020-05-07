#include "gva-utils.h"
#include <array>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GaussHermite {
namespace GVA {

template <class Type, class Fam>
integral_atomic<Type, Fam>&
integral_atomic<Type, Fam>::get_cached(unsigned const n){
  using output_T = integral_atomic<Type, Fam>;

  constexpr std::size_t const n_cache = GaussHermiteDataCachedMaxArg();
  if(n > n_cache or n == 0l)
    throw std::invalid_argument(
        "integral_atomic<Type, Fam>: invalid n (too large or zero)");

  static std::array<std::unique_ptr<output_T>, n_cache> cached_values;

  unsigned const idx = n - 1L;
  bool has_value = cached_values[idx].get();

  if(has_value)
    return *cached_values[idx];

#ifdef _OPENMP
#pragma omp critical (gvaCache)
{
#endif
  has_value = cached_values[idx].get();
  if(!has_value){
    std::unique_ptr<output_T> new_ptr(
        new output_T("integral_atomic<Type, Fam>", n));
    std::swap(cached_values[idx], new_ptr);
  }
#ifdef _OPENMP
}
#endif

  return *cached_values[idx];
}

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;

template class integral_atomic<double, mlogit_fam>;
template class integral_atomic<ADd   , mlogit_fam>;
template class integral_atomic<ADdd  , mlogit_fam>;
template class integral_atomic<ADddd , mlogit_fam>;

template class integral_atomic<double, probit_fam>;
template class integral_atomic<ADd   , probit_fam>;
template class integral_atomic<ADdd  , probit_fam>;
template class integral_atomic<ADddd , probit_fam>;

} // namespace GaussHermite
} // namespace GVA
