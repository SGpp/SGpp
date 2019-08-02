// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/MutexType.hpp>

namespace sgpp {
namespace base {

/**
 * Wrapper around MutexType which locks and unlocks upon
 * construction/destruction.
 *
 * Adopted from http://bisqwit.iki.fi/story/howto/openmp/#Locks.
 */
class ScopedLock {
 public:
  /**
   * Constructor, locks the MutexType object.
   *
   * @param m     MutexType object to be wrapped
   */
  explicit ScopedLock(MutexType& m) : mut(m), locked(true) { mut.lock(); }

  /**
   * Destructor, unlocks the MutexType object.
   */
  ~ScopedLock() { unlock(); }

  /**
   * Unlocks the MutexType object, if locked.
   */
  void unlock() {
    if (locked) {
      locked = false;
      mut.unlock();
    }
  }

  /**
   * Re-locks the MutexType object, if unlocked.
   */
  void lockAgain() {
    if (!locked) {
      mut.lock();
      locked = true;
    }
  }

 protected:
  /// underlying MutexType object
  MutexType& mut;
  /// whether the MutexType object is locked or not
  bool locked;

 private:
  /**
   * Custom copy constructor to prevent copying the lock.
   *
   * @param other     object to be copied
   */
  ScopedLock(const ScopedLock& other);

  /**
   * Custom assignment operator to prevent copying the lock.
   *
   * @param other     object to be assigned to
   */
  void operator=(const ScopedLock& other);
};
}  // namespace base
}  // namespace sgpp
