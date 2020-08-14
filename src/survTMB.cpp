#include "laplace.h"
#include "gva.h"
#include "snva.h"
#include "utils.h"

#include <cmath>
#include <algorithm>
#include <future>
#include <limits>

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace survTMB;

template<class Type>
Type objective_function<Type>::operator() ()
{
  SETUP_DATA;
  SETUP_DATA_CHECK;

#ifdef _OPENMP
  if(omp_get_max_threads() != n_threads)
    error("omp_get_max_threads() != n_threads");
#endif

  parallel_accumulator<Type> result(this);

  /* perform approximations using method descibed at
   *   https://github.com/kaskr/adcomp/issues/233#issuecomment-306032192
   *
   * each branch may contain further parameters or data objects */
  if(app_type == "Laplace"){
    PARAMETER_MATRIX(u);
    laplace(COMMON_CALL, u);
    return result;

  } else if(app_type =="GVA"){
    PARAMETER_VECTOR(theta_VA);
    DATA_INTEGER(n_nodes);
    GVA(COMMON_CALL, theta_VA, n_nodes);
    return result;

  }  else if(app_type == "SNVA"){
    PARAMETER_VECTOR(theta_VA);
    DATA_INTEGER(n_nodes);
    DATA_STRING(param_type)
    SNVA(COMMON_CALL, theta_VA, n_nodes, param_type);
    return result;

  }

  error("approximation method '%s' is not implemented", app_type.c_str());
  return Type(0);
}
