#ifndef SURVTMB_UTILS_H
#define SURVTMB_UTILS_H

#include "tmb_includes.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include <memory>

namespace survTMB {

inline unsigned get_rng_dim(unsigned const n_vcov_params){
  double const n(n_vcov_params);
  return std::lround(std::sqrt(.25 + 2 * n) - .5);
}

template<class Type>
unsigned get_rng_dim(vector<Type> x){
  return get_rng_dim(x.size());
}

/* Returns a covariance matix given a vector containing log standard
 * deviations, s, and lower triangular matrix elements theta. E.g. in 3x3

 \begin{align*}
 (\log \vec\sigma)^\top &= (s_1, s_2, s_3)
 & L &= \begin{pmatrix}
 1 & 0 & 0 \       \
 \theta_1 & 1 & 0 \\
 \theta_2 & \theta_3 & 1
 \end{pmatrix} \\
 \Sigma &= \text{diag}(\sigma)LL^\top\text{diag}(\sigma)
 \end{align*}
*/
template<class Type>
matrix<Type>
get_vcov_from_trian(Type const *vals, unsigned const dim){
  matrix<Type> out(dim, dim);
  out.setZero();

  /* fill the diagonal */
  for(unsigned i = 0; i < dim; ++i)
    out(i, i) = exp(*(vals + i));

  /* fill the off-diagonal */
  Type const * t = vals + dim;
  for(unsigned cl = 0; cl < dim; cl++)
    for(unsigned rw = cl + 1L; rw < dim; rw++)
      out(rw, cl) = out(rw, rw) * *t++;

  return out * out.transpose();
}

template<class Type>
matrix<Type>
get_vcov_from_trian(vector<Type> const &theta){
  unsigned const dim = get_rng_dim(theta);
  return get_vcov_from_trian(&theta[0], dim);
}

/* unlike objective_function::parallel_region this does not increase the
 * counter
 */
template<class Type>
bool is_my_region(objective_function<Type> const &o){
  if(o.current_parallel_region < 0 || o.selected_parallel_region < 0)
    /* Serial mode */
    return true;
  if(o.selected_parallel_region == o.current_parallel_region &&
     (!o.parallel_ignore_statements))
    return true;
  return false;
}

/* mock class to use instead when not using the TMB framework */
struct objective_mock {
  int selected_parallel_region = 0;

  inline bool is_my_region() const {
#ifdef _OPENMP
    return omp_get_thread_num() == selected_parallel_region;
#else
    return true;
#endif
  }

  inline bool parallel_region(){
    bool const ans = is_my_region();
#ifdef _OPENMP
    ++selected_parallel_region;
    if(selected_parallel_region >= omp_get_num_threads())
      selected_parallel_region = 0L;
#endif
    return ans;
  }
};

template<class Type>
class accumulator_mock {
  Type result = Type(0.);

public:
  std::unique_ptr<objective_mock> const obj;

  accumulator_mock(): obj(new objective_mock()) { }

  inline void operator+=(Type x){
    if(obj->parallel_region())
      result += x;
  }
  inline void operator-=(Type x){
    if(obj->parallel_region())
      result -= x;
  }
  operator Type(){
    return result;
  }
};

inline bool is_my_region(objective_mock const &o){
  return o.is_my_region();
}

} // namespace survTMB

#endif
