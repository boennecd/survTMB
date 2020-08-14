#ifndef CLEAR_MEM_H
#define CLEAR_MEM_H
#include "tmb_includes.h"

template<class Type>
void track_atomic(CppAD::atomic_base<Type>*);

template<template<class> class Cont, class Type>
Cont<Type>* track_atomic(Cont<Type> *ptr){
  track_atomic((CppAD::atomic_base<Type>*)ptr);
  return ptr;
}

template<class Type>
void clear_atomics();

class clearable {
public:
  virtual ~clearable() = default;
  virtual void clear() { };
};

void add_clearable(clearable*);
void clear_clearables();

#endif
