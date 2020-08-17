#ifndef SURVTMB_UTILS_H
#define SURVTMB_UTILS_H

#include "tmb_includes.h"
#include "fast-commutation.h"
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

/* Returns the covariance matrix from a log Cholesky decomposition
*/
template<class Type>
class get_vcov_from_trian_atomic : public CppAD::atomic_base<Type> {
public:
  get_vcov_from_trian_atomic(char const *name):
  CppAD::atomic_base<Type>(name) {
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }

  template<class T>
  static void comp(T const *x, T * y, size_t const n_arg){
    // TODO: can be done smarter...
    typename
    Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >
      out(y, n_arg, n_arg);
    out.setZero();

    for(size_t c = 0; c < n_arg; ++c){
      for(size_t r = c; r < n_arg; ++r)
        out(r, c) = *x++;
    }

    // out * out.t()
    for(size_t c = n_arg; c-- > 0;){
      for(size_t r = c + 1L; r-- > 0;){
        Type new_term(0);
        for(size_t k = 0; k <= r; ++k)
          new_term += out(c, k) * out(r, k);
        out(c, r) = new_term;
        out(r, c) = new_term;
      }
    }
  }

  virtual bool forward(std::size_t p, std::size_t q,
                       const CppAD::vector<bool> &vx,
                       CppAD::vector<bool> &vy,
                       const CppAD::vector<Type> &tx,
                       CppAD::vector<Type> &ty){
    if(q > 0)
      return false;

    size_t const n = std::lround(
      .5 * (std::sqrt(8. * static_cast<double>(tx.size()) + 1.) - 1.));
    comp(&tx[0], &ty[0], n);

    /* set variable flags */
    if (vx.size() > 0) {
      bool anyvx = false;
      for (std::size_t i = 0; i < vx.size() and !anyvx; i++)
        anyvx |= vx[i];
      for (std::size_t i = 0; i < vy.size(); i++)
        vy[i] = anyvx;
    }

    return true;
  }

  virtual bool reverse(std::size_t q, const CppAD::vector<Type> &tx,
                       const CppAD::vector<Type> &ty,
                       CppAD::vector<Type> &px,
                       const CppAD::vector<Type> &py){
    if(q > 0)
      return false;

    size_t const n = std::lround(
      .5 * (std::sqrt(8. * static_cast<double>(tx.size()) + 1.) - 1.)),
                nn = n * n;

    size_t const * const com_vec = get_commutation_unequal_vec_cached(n);
    size_t bri(0L),
           brl(0L);

    for(size_t i = 0; i < px.size(); ++i)
      px[i] = Type(0.);

    for(size_t i = 0L; i < nn;
        bri = bri + 1L >= n ? 0L : bri + 1L,
        brl = bri == 0L ? brl + 1L : brl,
        ++i){
      size_t const idx_l = *(com_vec + i);
      Type const lfs = py[idx_l] + py[i];

      size_t idx_k(bri);
      for(size_t k = 0; k <= bri; idx_k += n - k, ++k){
        if(brl >= k){
          Type const rfs = tx[
            brl - k /* row index */ + idx_k - bri /* col index */];
          px[idx_k - k] += lfs * rfs;
        }
      }
    }

    return true;
  }

  virtual bool rev_sparse_jac(size_t q, const CppAD::vector<bool>& rt,
                              CppAD::vector<bool>& st) {
    bool anyrt = false;
    for (std::size_t i = 0; i < rt.size() and !anyrt; i++)
      anyrt |= rt[i];
    for (std::size_t i = 0; i < st.size(); i++)
      st[i] = anyrt;
    return true;
  }

  static get_vcov_from_trian_atomic<Type>& get_cached();
};

template<class Type>
matrix<AD<Type> >
get_vcov_from_trian(AD<Type> const *vals, unsigned const dim){
  size_t const n_ele = (dim * (dim + 1L)) / 2L;
  matrix<AD<Type> > out(dim, dim);
  CppAD::vector<AD<Type> > x_vals(n_ele);

  for(size_t i = 0; i < n_ele; ++i)
    x_vals[i] = *(vals + i);
  size_t cc = 0L;
  for(size_t i = dim + 1L; i-- > 1; cc += i)
    x_vals[cc] = exp(x_vals[cc]);

  static get_vcov_from_trian_atomic<Type> &functor =
    get_vcov_from_trian_atomic<Type>::get_cached();

  typename Eigen::Map<Eigen::Matrix<AD<Type> ,Eigen::Dynamic,1> >
    tx(x_vals.data(), n_ele),
    ty(out.data(), dim * dim);

  functor(tx, ty);

  return out;
}

inline matrix<double>
get_vcov_from_trian(double const *vals, unsigned const dim){
  matrix<double> out(dim, dim);
  size_t const n_ele = (dim * (dim + 1L)) / 2L;
  vector<double> x_vals(n_ele);
  for(size_t i = 0; i < n_ele; ++i)
    x_vals[i] = *(vals + i);
  size_t cc = 0L;
  for(size_t i = dim + 1L; i-- > 1; cc += i)
    x_vals[cc] = exp(x_vals[cc]);

  double *p_out = out.data();
  get_vcov_from_trian_atomic<double>::comp(x_vals.data(), p_out, dim);

  return out;
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
class objective_mock {
#ifdef _OPENMP
  int const my_num = omp_get_thread_num(),
         n_regions = omp_get_num_threads();
#endif

public:
  int selected_parallel_region = 0;

  inline bool is_my_region() const {
#ifdef _OPENMP
    return my_num == selected_parallel_region;
#else
    return true;
#endif
  }

  inline bool parallel_region(){
    bool const ans = is_my_region();
#ifdef _OPENMP
    ++selected_parallel_region;
    if(selected_parallel_region >= n_regions)
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
