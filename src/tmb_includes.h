#ifndef TMB_INCLUDES_H
#define TMB_INCLUDES_H

#ifdef INCLUDE_RCPP
/* Rcpp */
#include <RcppArmadillo.h>

/* Rcpp defines `R_NO_REMAP` which causes errors when using TMB. Thus, we
* define a few macros. */
#ifndef findVar
#define findVar	Rf_findVar
#endif
#ifndef install
#define install	Rf_install
#endif
#ifndef error
#define error	Rf_error
#endif

/* TMB redefines a few macro. Thus, we undefine them first */
#ifdef R_D_val
#undef R_D_val
#endif
#ifdef R_D_qIv
#undef R_D_qIv
#endif
#ifdef R_D_exp
#undef R_D_exp
#endif
#ifdef R_D_log
#undef R_D_log
#endif
#ifdef R_D_Clog
#undef R_D_Clog
#endif
#ifdef R_D_LExp
#undef R_D_LExp
#endif
#ifdef R_DT_qIv
#undef R_DT_qIv
#endif
#ifdef R_DT_CIv
#undef R_DT_CIv
#endif
#ifdef R_Q_P01_check
#undef R_Q_P01_check
#endif
#ifdef R_Q_P01_boundaries
#undef R_Q_P01_boundaries
#endif
#ifdef R_D_negInonint
#undef R_D_negInonint
#endif
#ifdef R_D_nonint_check
#undef R_D_nonint_check
#endif
#endif // ifdef WITH_RCPP

/* TMB */
#define WITH_LIBTMB

#ifdef TMB_PRECOMPILE
#undef TMB_PRECOMPILE
#endif

#ifdef CSKIP
#undef CSKIP
#endif

#include <TMB.hpp>

namespace CppAD {
inline double Value(double const x){
  return x;
}
}

template<class Type>
Type vec_dot(vector<Type> const &x, vector<Type> const &y){
#ifdef DO_CHECKS
  if(x.size() != y.size())
    throw std::invalid_argument("vec_dot<Type>: dimension do not match");
#endif
  Type out(0.);
  for(int i = 0; i < x.size(); ++i)
    out += x[i] * y[i];
  return out;
}

template<class Type>
Type vec_dot
  (Eigen::Map<const Eigen::Matrix<Type,Eigen::Dynamic,1> > &x,
   vector<Type> const &y){
#ifdef DO_CHECKS
  if(x.size() != y.size())
    throw std::invalid_argument("vec_dot<Type>: dimension do not match");
#endif
  Type out(0.);
  for(int i = 0; i < x.size(); ++i)
    out += x[i] * y[i];
  return out;
}

#ifdef INCLUDE_RCPP
inline double vec_dot(vector<double> const &x, arma::vec const &y){
#ifdef DO_CHECKS
  if(x.size() != y.n_elem)
    throw std::invalid_argument("vec_dot: dimension do not match");
#endif
  double out(0.);
  for(int i = 0; i < x.size(); ++i)
    out += x[i] * y[i];
  return out;
}
#endif

template<class Type>
Type quad_form
  (vector<Type> const &x, matrix<Type> const &A, vector<Type> const &y){
#ifdef DO_CHECKS
  if(x.size() != A.rows() or y.size() != A.cols())
    throw std::invalid_argument("quad_form<Type>: dimension do not match");
#endif
  Type out(0.);
  for(int j = 0; j < A.cols(); ++j)
    for(int i = 0; i < A.rows(); ++i)
      out += A(i, j) * x[i] * y[j];

  return out;
}

template<class Type>
Type quad_form_sym(vector<Type> const &x, matrix<Type> const &A){
#ifdef DO_CHECKS
  if(x.size() != A.rows() or x.size() != A.cols())
    throw std::invalid_argument(
        "quad_form_sym<Type>: dimension do not match");
#endif
  Type out(0.);
  Type const TWO(2.);
  for(int j = 0; j < A.cols(); ++j){
    for(int i = 0; i < j; ++i)
      out += TWO * A(i, j) * x[i] * x[j];
    out += A(j, j) * x[j] * x[j];
  }

  return out;
}

template<class Type>
Type mat_mult_trace(matrix<Type> const &X, matrix<Type> const &Y){
#ifdef DO_CHECKS
  if(X.rows() != Y.cols() or X.cols() != Y.rows())
    throw std::invalid_argument(
        "mat_mult_trace<Type>: dimension do not match");
#endif
  Type out(0.);
  for(int j = 0; j < X.cols(); ++j)
    for(int i = 0; i < X.rows(); ++i)
      out += X(i, j) * Y(i, j);

  return out;
}

#endif
