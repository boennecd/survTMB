#ifndef FAST_COMMUTATION_H
#define FAST_COMMUTATION_H
#include <memory>

std::unique_ptr<size_t[]> get_commutation_unequal_vec
  (unsigned const, unsigned const, bool const);

size_t const * get_commutation_unequal_vec_cached(unsigned const);

#endif
