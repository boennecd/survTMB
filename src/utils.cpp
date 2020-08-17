#include "utils.h"
#include <array>
#include <memory.h>
#include <string.h>
#include "clear-mem.h"

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
#ifdef _OPENMP
#pragma omp master
  track_atomic(&get_vcov_from_trian_d);
#else
{
  static bool has_added = false;
  if(!has_added){
    track_atomic(&get_vcov_from_trian_d);
    has_added = true;
  }
}
#endif

  return get_vcov_from_trian_d;
}


template<>
get_vcov_from_trian_atomic<ADd>&
get_vcov_from_trian_atomic<ADd>::get_cached(){
#ifdef _OPENMP
#pragma omp master
  track_atomic(&get_vcov_from_trian_ADd);
#else
{
  static bool has_added = false;
  if(!has_added){
    track_atomic(&get_vcov_from_trian_ADd);
    has_added = true;
  }
}
#endif

  return get_vcov_from_trian_ADd;
}

template<>
get_vcov_from_trian_atomic<ADdd>&
get_vcov_from_trian_atomic<ADdd>::get_cached(){
#ifdef _OPENMP
#pragma omp master
  track_atomic(&get_vcov_from_trian_ADdd);
#else
{
  static bool has_added = false;
  if(!has_added){
    track_atomic(&get_vcov_from_trian_ADdd);
    has_added = true;
  }
}
#endif

  return get_vcov_from_trian_ADdd;
}

template<>
get_vcov_from_trian_atomic<ADddd>&
get_vcov_from_trian_atomic<ADddd>::get_cached(){
#ifdef _OPENMP
#pragma omp master
  track_atomic(&get_vcov_from_trian_ADddd);
#else
{
  static bool has_added = false;
  if(!has_added){
    track_atomic(&get_vcov_from_trian_ADddd);
    has_added = true;
  }
}
#endif

  return get_vcov_from_trian_ADddd;
}

} // namespace

