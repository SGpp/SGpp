// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_IDENTITYPOINTORDERING_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_IDENTITYPOINTORDERING_HPP_

#include <sgpp/combigrid/grid/growth/AbstractGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/ordering/AbstractPointOrdering.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Point ordering class that uses the identity mapping to compute a level-independent index (i.e. it
 * can be used for non-nested settings or for points which are already ordered correctly, e.g. Leja
 * points). The growth strategy can be configured via a constructor parameter. If isSorted is set to
 * false in the constructor, the sorted permutation is determined using a sorting algorithm.
 */
class IdentityPointOrdering : public AbstractPointOrdering {
  std::shared_ptr<AbstractGrowthStrategy> growthStrategy;
  bool isSorted;

 public:
  IdentityPointOrdering(std::shared_ptr<AbstractGrowthStrategy> growthStrategy,
                        bool isSorted = false);

  virtual ~IdentityPointOrdering();

  virtual size_t convertIndex(size_t level, size_t numPoints, size_t index);

  virtual size_t numPoints(size_t level);

  virtual std::shared_ptr<AbstractPermutationIterator> getSortedPermutationIterator(
      size_t level, std::vector<double> const &points, size_t numPoints);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_IDENTITYPOINTORDERING_HPP_ */
