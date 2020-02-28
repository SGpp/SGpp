// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

namespace sgpp {
namespace base {

#ifdef _OPENMP

/**
 * Wrapper for OpenMP nested locks.
 * Once locked, other threads block when trying to lock the same MutexType.
 * However, the MutexType can be locked multiple times by the same
 * thread without blocking.
 * In this case, it has to be unlocked the same number of times to
 * let the other threads continuing execution.
 *
 * Adopted from http://bisqwit.iki.fi/story/howto/openmp/#Locks.
 */
class MutexType {
 public:
  /**
   * Constructor.
   */
  MutexType() { omp_init_nest_lock(&_lock); }

  /**
   * Custom copy constructor to prevent copying the lock.
   *
   * @param other     object to be copied
   */
  MutexType(const MutexType& other) {
    (void)other;
    omp_init_nest_lock(&_lock);
  }

  /**
   * Custom assignment operator to prevent copying the lock
   *
   * @param other     object to be assigned to
   * @return          *this
   */
  MutexType& operator=(const MutexType& other) {
    (void)other;
    return *this;
  }

  /**
   * Destructor.
   */
  ~MutexType() { omp_destroy_nest_lock(&_lock); }

  /**
   * Lock the MutexType.
   */
  void lock() { omp_set_nest_lock(&_lock); }

  /**
   * Unlock the MutexType.
   */
  void unlock() { omp_unset_nest_lock(&_lock); }

 protected:
  /// OpenMP lock
  omp_nest_lock_t _lock;
};

#else

// dummy definition in case SG++ is compiled without OpenMP
class MutexType {
 public:
  void lock() {}

  void unlock() {}
};

#endif /* _OPENMP */
}  // namespace base
}  // namespace sgpp
