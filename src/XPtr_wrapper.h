#ifndef XPTR_WRAPPER_H
#define XPTR_WRAPPER_H
#define INCLUDE_RCPP
#include "clear-mem.h"

template<typename T>
class XPtr_wrapper final : public clearable {
  Rcpp::XPtr<T> ptr;

public:
  XPtr_wrapper(T * obj): ptr(obj, false) { }

  void clear() {
    ptr.release();
  }

  explicit operator Rcpp::XPtr<T>() const {
    return ptr;
  }
};

#endif
