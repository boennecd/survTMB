#include "tmb_includes.h"
#include "clear-mem.h"

template<class Type>
class vec_dot_atomic : public CppAD::atomic_base<Type> {
public:
  vec_dot_atomic(char const *name):
  CppAD::atomic_base<Type>(name) {
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }

  template<class T>
  static T comp(T const *x, T const *y, size_t const n){
    Type out(0.);
    for(size_t i = 0; i < n; ++i, ++x, ++y)
      out += *x * *y;

    return out;
  }

  virtual bool forward(std::size_t p, std::size_t q,
                       const CppAD::vector<bool> &vx,
                       CppAD::vector<bool> &vy,
                       const CppAD::vector<Type> &tx,
                       CppAD::vector<Type> &ty){
    if(q > 0)
      return false;

    size_t const n = tx.size() / 2L;
    Type const * xx = &tx[0],
               * yy = xx + n;
    ty[0] = comp(xx, yy, n);

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

    size_t const n = tx.size() / 2L;
    Type const * xx = &tx[0],
               * yy = xx + n;

    size_t j = n;
    for(size_t i = 0L; i < n; ++i, ++j, ++xx, ++yy){
      px[i] = *yy * py[0];
      px[j] = *xx * py[0];
    }

    return true;
  }

  virtual bool rev_sparse_jac(size_t q, const CppAD::vector<bool>& rt,
                              CppAD::vector<bool>& st) {
    bool anyrt = false;
    for (std::size_t i = 0; i < rt.size(); i++)
      anyrt |= rt[i];
    for (std::size_t i = 0; i < st.size(); i++)
      st[i] = anyrt;
    return true;
  }
};

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;

static vec_dot_atomic<double> vec_dot_d     = vec_dot_atomic<double>("vec_dot_atomic<double>");
static vec_dot_atomic<ADd   > vec_dot_ADd   = vec_dot_atomic<ADd   >("vec_dot_atomic<AD<double> >");
static vec_dot_atomic<ADdd  > vec_dot_ADdd  = vec_dot_atomic<ADdd  >("vec_dot_atomic<AD<AD<double> > >");

template<>
double vec_dot(double const *x, double const *y, size_t const n){
  return vec_dot_atomic<double>::comp(x, y, n);
}

#define VEC_DOT_SPEC(OBJ_NAME, TNAME)                          \
template<>                                                     \
TNAME vec_dot(TNAME const *x, TNAME const *y, size_t const n){ \
  CppAD::vector<TNAME> xx(2L * n),                             \
                       yy(1L);                                 \
                                                               \
  size_t i = 0;                                                \
  for(; i < n; ++i, ++x)                                       \
    xx[i] = *x;                                                \
  size_t const e2 = 2L * n;                                    \
  for(; i < e2; ++i, ++y)                                      \
    xx[i] = *y;                                                \
                                                               \
  OBJ_NAME(xx, yy, n);                                         \
                                                               \
  return yy[0];                                                \
}

VEC_DOT_SPEC(vec_dot_d   , ADd)
VEC_DOT_SPEC(vec_dot_ADd , ADdd)
VEC_DOT_SPEC(vec_dot_ADdd, ADddd)

#undef VEC_DOT_SPEC

template<class Type>
class quad_form_atomic : public CppAD::atomic_base<Type> {
public:
  quad_form_atomic(char const *name):
  CppAD::atomic_base<Type>(name) {
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);
  }

  template<class T>
  static T comp(T const *x, T const *a, T const *y, size_t const n){
    typename
    Eigen::Map<const Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic> >
      A(a, n, n);

    Type out(0.);
    for(size_t j = 0; j < n; ++j)
      for(size_t i = 0; i < n; ++i)
        out += A(i, j) * *(x + i) * *(y + j);

    return out;
  }

  virtual bool forward(std::size_t p, std::size_t q,
                       const CppAD::vector<bool> &vx,
                       CppAD::vector<bool> &vy,
                       const CppAD::vector<Type> &tx,
                       CppAD::vector<Type> &ty){
    if(q > 0)
      return false;

    size_t const n =
      std::lround(std::sqrt(static_cast<double>(tx.size()) + 1.) - 1.);
    Type const * x = &tx[0],
               * a = x + n,
               * y = a + n * n;
    ty[0] = comp(x, a, y, n);

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

    size_t const n =
      std::lround(std::sqrt(static_cast<double>(tx.size()) + 1.) - 1.);
    Type const * x = &tx[0],
               * a = x + n,
               * y = a + n * n;
    Type * dx = &px[0],
         * da = dx + n,
         * dy = da + n * n;

    typename
    Eigen::Map<const Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic> >
      A(a, n, n);

    for(size_t i = 0; i < px.size(); ++i)
      px[i] = Type(0.);

    for(size_t j = 0; j < n; ++j){
      size_t const jn = j * n;
      for(size_t i = 0; i < n; ++i){
        *(dx + i)      += A(i, j) *  *(y + j);
        // TODO: there is an assumption here about column major order
        *(da + i + jn) += *(x + i) * *(y + j);
        *(dy + j)      += A(i, j) * *(x + i);
      }
    }

    for(size_t i = 0; i < px.size(); ++i)
      px[i] *= py[0];

    return true;
  }

  virtual bool rev_sparse_jac(size_t q, const CppAD::vector<bool>& rt,
                              CppAD::vector<bool>& st) {
    bool anyrt = false;
    for (std::size_t i = 0; i < rt.size(); i++)
      anyrt |= rt[i];
    for (std::size_t i = 0; i < st.size(); i++)
      st[i] = anyrt;
    return true;
  }
};

static quad_form_atomic<double> quad_form_d     = quad_form_atomic<double>("quad_form_atomic<double>");
static quad_form_atomic<ADd   > quad_form_ADd   = quad_form_atomic<ADd   >("quad_form_atomic<AD<double> >");
static quad_form_atomic<ADdd  > quad_form_ADdd  = quad_form_atomic<ADdd  >("quad_form_atomic<AD<AD<double> > >");

template<>
double quad_form(double const *x, double const *a, double const *y,
                 size_t const n){
  return quad_form_atomic<double>::comp(x, a, y, n);
}

#define QUAD_FORM_SPEC(OBJ_NAME, TNAME)                        \
template<>                                                     \
TNAME quad_form(TNAME const *x, TNAME const *a, TNAME const *y,\
                size_t const n){                               \
  CppAD::vector<TNAME> arg(2L * n + n * n),                    \
                       out(1L);                                \
                                                               \
  size_t i = 0L;                                               \
  for(; i < n; ++i, ++x)                                       \
    arg[i] = *x;                                               \
  size_t const e1 = n + n * n;                                 \
  for(; i < e1; ++i, ++a)                                      \
    arg[i] = *a;                                               \
  size_t e2 = 2L * n + n * n;                                  \
  for(; i < e2; ++i, ++y)                                      \
    arg[i] = *y;                                               \
                                                               \
  OBJ_NAME(arg, out);                                          \
                                                               \
  return out[0L];                                              \
}

QUAD_FORM_SPEC(quad_form_d   , ADd)
QUAD_FORM_SPEC(quad_form_ADd , ADdd)
QUAD_FORM_SPEC(quad_form_ADdd, ADddd)


void add_atomics_to_be_cleared(){
  static bool have_been_added = false;
  if(have_been_added)
    return;

  track_atomic(&vec_dot_d);
  track_atomic(&vec_dot_ADd);
  track_atomic(&vec_dot_ADdd);

  track_atomic(&quad_form_d);
  track_atomic(&quad_form_ADd);
  track_atomic(&quad_form_ADdd);
  have_been_added = true;
}

#undef QUAD_FORM_SPEC
