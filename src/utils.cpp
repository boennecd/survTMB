#include "utils.h"
#include <array>
#include <memory.h>
#include <string.h>

namespace survTMB {

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;

static get_vcov_from_trian_atomic<double> get_vcov_from_trian_d
  ("get_vcov_from_trian_atomic<double>");
static get_vcov_from_trian_atomic<ADd   > get_vcov_from_trian_ADd
  ("get_vcov_from_trian_atomic<ADd>");
static get_vcov_from_trian_atomic<ADdd  > get_vcov_from_trian_ADdd
  ("get_vcov_from_trian_atomic<ADdd>");
static get_vcov_from_trian_atomic<ADddd > get_vcov_from_trian_ADddd
  ("get_vcov_from_trian_atomic<ADddd>");

template<>
get_vcov_from_trian_atomic<double>&
get_vcov_from_trian_atomic<double>::get_cached(){
  return get_vcov_from_trian_d;
}


template<>
get_vcov_from_trian_atomic<ADd>&
get_vcov_from_trian_atomic<ADd>::get_cached(){
  return get_vcov_from_trian_ADd;
}

template<>
get_vcov_from_trian_atomic<ADdd>&
get_vcov_from_trian_atomic<ADdd>::get_cached(){
  return get_vcov_from_trian_ADdd;
}

template<>
get_vcov_from_trian_atomic<ADddd>&
get_vcov_from_trian_atomic<ADddd>::get_cached(){
  return get_vcov_from_trian_ADddd;
}

} // namespace

