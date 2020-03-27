#include "laplace.h"
#include "gva.h"
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
  /* common data objects and parameters */
  DATA_STRING(app_type);
  DATA_VECTOR(tobs);

  DATA_VECTOR(event);
  DATA_MATRIX(X);
  DATA_MATRIX(XD);
  DATA_MATRIX(Z);
  DATA_IVECTOR(grp);
  DATA_STRING(link);
  DATA_IVECTOR(grp_size);
  DATA_INTEGER(n_threads)

  /* They are marked as parameter such that the user can change them
   * later */
  PARAMETER(eps);
  PARAMETER(kappa);

  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(theta);

  /* checks */
  {
    unsigned const n = tobs.size();
    auto check_rows = [n](matrix<Type> const &x, char const *msg){
      if(x.rows() != n)
        error(msg);
    };
    check_rows(X , "invalid 'X'");
    check_rows(XD, "invalid 'XD'");
    check_rows(Z , "invalid 'Z'");
    if(n != event.size())
      error("invalid 'event'");
    if(n != grp.size())
      error("invalid 'grp'");

    if(b.size()  != X.cols())
      error("invalid 'b'");
    if(XD.cols() != X.cols())
      error("invalid 'XD' (# columns)");

    if(Z.cols() != get_rng_dim(theta))
      error("invalid 'Z' (# columns: %d %d)", Z.cols(), get_rng_dim(theta));

    std::size_t grp_size_sum = 0;
    for(int i = 0; i < grp_size.size(); ++i)
      grp_size_sum += grp_size[i];
    if(n != grp_size_sum)
      error("invalid 'grp_size'");

    int g_max = 0;
    for(int i = 0; i < grp.size(); ++i)
      g_max = std::max(grp[i], g_max);
    if(g_max + 1L != grp_size.size())
      error("invalid 'grp_size' (does not match with number of groups)");

    for(int i = 1L; i < grp.size(); ++i){
      if(grp[i - 1L] > grp[i])
        error("'grp' is not sorted");
      if(grp[i - 1L] - grp[i] < -1L)
        error("too big gap in 'grp'");
    }
  }

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
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

  }

  error("approximation method '%s' is not implemented", app_type.c_str());
  return Type(0);
}
