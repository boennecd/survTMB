#ifndef CONVERT_EIGEN_ARMA_H
#define CONVERT_EIGEN_ARMA_H
#define INCLUDE_RCPP
#include "tmb_includes.h"

inline arma::vec
vec_eigen_arma(vector<double> const &ev){
  arma::vec out(ev.size());
  for(unsigned i = 0; i < ev.size(); ++i)
    out[i] = ev[i];
  return out;
}

template<typename Type = double>
vector<Type>
vec_eigen_arma(arma::vec const &av){
  vector<Type> out(av.n_elem);
  for(unsigned i = 0; i < av.n_elem; ++i)
    out[i] = Type(av[i]);
  return out;
}

template<typename Type = double>
matrix<Type>
mat_eigen_arma(arma::mat const &am){
  matrix<Type> out(am.n_rows, am.n_cols);
  for(unsigned j = 0; j < am.n_cols; j++)
    for(unsigned i = 0; i < am.n_rows; ++i)
      out(i, j) = Type(am(i, j));
  return out;
}

#endif
