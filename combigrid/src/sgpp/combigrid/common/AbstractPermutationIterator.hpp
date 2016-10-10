// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_

#include <sgpp/globaldef.hpp>

#include <cstddef>

#include <memory>

namespace sgpp {
namespace combigrid {

/**
 * Iterates over the elements of a permutation.
 * It is assumed that the user already knows the number of elements of the permutation.
 * This class is used internally to obtain grid points in a sorted order, i.e. from left to right,
 * by providing the indices in the correct order.
 */
class AbstractPermutationIterator {
 public:
  virtual ~AbstractPermutationIterator();

  /**
   * Sets the iterator back to the start
   */
  virtual void reset() = 0;

  /**
   * Current value of the permutation.
   */
  virtual size_t value() = 0;

  /**
   * Go to next element of the permutation.
   */
  virtual void moveToNext() = 0;

  virtual std::shared_ptr<AbstractPermutationIterator> clone() = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_ */
