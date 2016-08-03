// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_SORTEDPERMUTATIONITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_SORTEDPERMUTATIONITERATOR_HPP_

#include <sgpp/combigrid/common/AbstractPermutationIterator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class SortedPermutationIterator : public AbstractPermutationIterator {
  std::vector<size_t> permutation;
  size_t currentIndex;

  SortedPermutationIterator(std::vector<size_t> const &permutation, size_t currentIndex);

 public:
  SortedPermutationIterator(std::vector<double> const &points, size_t numPoints);
  virtual ~SortedPermutationIterator();

  /**
   * Sets the iterator back to the start
   */
  virtual void reset();

  virtual size_t value();

  virtual void moveToNext();

  virtual std::shared_ptr<AbstractPermutationIterator> clone();
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_SORTEDPERMUTATIONITERATOR_HPP_ */
