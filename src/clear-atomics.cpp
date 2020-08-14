#include "clear-mem.h"
#include <vector>
#include "memory.h"

using ADd   = CppAD::AD<double>;
using ADdd  = CppAD::AD<CppAD::AD<double> >;
using ADddd = CppAD::AD<CppAD::AD<CppAD::AD<double> > >;

namespace {
static std::vector<CppAD::atomic_base<double>* > vec_atomics_d;
static std::vector<CppAD::atomic_base<ADd   >* > vec_atomics_ADd;
static std::vector<CppAD::atomic_base<ADdd  >* > vec_atomics_ADdd;
static std::vector<CppAD::atomic_base<ADddd >* > vec_atomics_ADddd;

static std::vector<std::unique_ptr<clearable> > vec_clearable;
} // namespace

#define TRACK_ATOMIC_SPEC(TNAM, OBJ)                           \
template<>                                                     \
void track_atomic(CppAD::atomic_base<TNAM> *ptr){              \
  OBJ.emplace_back(ptr);                                       \
}

TRACK_ATOMIC_SPEC(double, vec_atomics_d)
TRACK_ATOMIC_SPEC(ADd   , vec_atomics_ADd)
TRACK_ATOMIC_SPEC(ADdd  , vec_atomics_ADdd)
TRACK_ATOMIC_SPEC(ADddd , vec_atomics_ADddd)

#undef TRACK_ATOMIC_SPEC

#define CLEAR_ATOMICS_SPEC(TNAM, OBJ)                          \
template<>                                                     \
void clear_atomics<TNAM>(){                                    \
  for(auto x : OBJ)                                            \
    if(x)                                                      \
      x->clear();                                              \
}

CLEAR_ATOMICS_SPEC(double, vec_atomics_d)
CLEAR_ATOMICS_SPEC(ADd   , vec_atomics_ADd)
CLEAR_ATOMICS_SPEC(ADdd  , vec_atomics_ADdd)
CLEAR_ATOMICS_SPEC(ADddd , vec_atomics_ADddd)

#undef CLEAR_ATOMICS_SPEC

void add_clearable(clearable *new_ele){
  vec_clearable.emplace_back(new_ele);
}

void clear_clearables(){
  for(auto &e : vec_clearable)
    if(e)
      e->clear();
  vec_clearable.clear();
}
