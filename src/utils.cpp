#include "utils.h"
#include <array>
#include <memory.h>
#include <string.h>

namespace survTMB {

template<class Type>
get_vcov_from_trian_atomic<Type>&
get_vcov_from_trian_atomic<Type>::get_cached(size_t const n_arg){
  using output_T = get_vcov_from_trian_atomic<Type>;

  constexpr std::size_t const n_cache = 1000L;
  if(n_arg > n_cache or n_arg == 0l)
    throw std::invalid_argument(
        "get_vcov_from_trian_atomic<Type>::get_cached: invalid n, " + std::to_string(n_arg) + ", (too large or zero)");

  static std::array<std::unique_ptr<output_T>, n_cache> cached_values;

  unsigned const idx = n_arg - 1L;
  bool has_value = cached_values[idx].get();

  if(has_value)
    return *cached_values[idx];

#ifdef _OPENMP
  if(CppAD::thread_alloc::in_parallel())
    throw std::runtime_error("get_vcov_from_trian_atomic<Type>::get_cached called in parallel mode");
#endif

  cached_values[idx].reset(new output_T("get_vcov_from_trian_atomic<Type>", n_arg));

  return *cached_values[idx];
}

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;

template class get_vcov_from_trian_atomic<double>;
template class get_vcov_from_trian_atomic<ADd   >;
template class get_vcov_from_trian_atomic<ADdd  >;
template class get_vcov_from_trian_atomic<ADddd >;

} // namespace

