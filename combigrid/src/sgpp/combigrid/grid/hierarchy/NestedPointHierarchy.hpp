// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_HIERARCHY_NESTEDPOINTHIERARCHY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_HIERARCHY_NESTEDPOINTHIERARCHY_HPP_

#include <sgpp/combigrid/common/AbstractPermutationIterator.hpp>
#include <sgpp/combigrid/grid/distribution/AbstractPointDistribution.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/AbstractPointOrdering.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * PointHierarchy class for nested point hierarchies (computes the nested points only once), see
 * also AbstractPointHierarchy
 */
class NestedPointHierarchy : public AbstractPointHierarchy {
  std::vector<size_t> numPointsPerLevel;
  std::vector<double> points;
  std::vector<std::shared_ptr<AbstractPermutationIterator>> permutationIterators;

  std::shared_ptr<AbstractPointDistribution> pointDistribution;
  std::shared_ptr<AbstractPointOrdering> pointOrdering;

 public:
  NestedPointHierarchy(std::shared_ptr<AbstractPointDistribution> pointDistribution,
                       std::shared_ptr<AbstractPointOrdering> pointOrdering);

  virtual ~NestedPointHierarchy();

  /**
   * @return Returns the grid point for the given level and index. (0 <= index <
   * getNumPoints(level)). In the case of NestedPointHierarchy, this also works if the argument
   * level is smaller than the real level. This is necessary because the CombigridTreeStorage loses
   * the level information on the way to calling this method.
   */
  virtual double getPoint(size_t level, size_t index);

  virtual std::vector<double> &computePoints(size_t level);
  virtual std::vector<double> getPoints(size_t level, bool sorted);

  virtual size_t getNumPoints(size_t level);

  virtual bool isNested();

  virtual std::shared_ptr<AbstractPermutationIterator> getSortedPermutationIterator(size_t level);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_HIERARCHY_NESTEDPOINTHIERARCHY_HPP_ */
