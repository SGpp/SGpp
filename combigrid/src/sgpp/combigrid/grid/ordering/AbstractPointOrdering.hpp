/*
 * AbstractPointOrdering.hpp
 *
 *  Created on: 18.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_ABSTRACTPOINTORDERING_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_ABSTRACTPOINTORDERING_HPP_

#include "../../common/AbstractPermutationIterator.hpp"
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <vector>
#include <memory>

namespace sgpp{
namespace combigrid {

/**
 * An AbstractPointOrdering may (via convertIndex()) define a reordering of the points so that the same points in different levels have the same indices.
 * Furthermore, it provides a SortedPermutationIterator to traverse the points in ascending order.
 *
 * Since the necessity of reordering depends on whether the points are nested and the nesting in turn may depend on the growth strategy,
 * an AbstractPointOrdering also includes a method to compute the number of points for a level.
 * This might be implemented using an AbstractGrowthStrategy-Object,
 * but can also be fixed because the ordering depends on the growth strategy.
 *
 * As an example, consider the UniformPointDistribution combined with an ExponentialLevelorderPointOrdering.
 * In Level 0 there is one point: 0.5
 * In Level 1 there are three points: 0.0, 0.5, 1.0
 * In order to exploit nesting, 0.5 needs to have the same index in all levels.
 * Thus, in level 1, the index mapping is 0 -> 1, 1 -> 0, 2 -> 2.
 * Since the original points were sorted but the remapped ones are not, the SortedPermutationIterator has to implement the inverse mapping.
 */
class AbstractPointOrdering {
public:
	virtual ~AbstractPointOrdering();

	/**
	 * @param level The level where the index is considered
	 * @param numPoints The number of points at the level
	 * @param index The index of the concrete point in the level
	 */
	virtual size_t convertIndex(size_t level, size_t numPoints, size_t index) = 0;

	virtual size_t numPoints(size_t level) = 0;

	virtual std::shared_ptr<AbstractPermutationIterator> getSortedPermutationIterator(size_t level, std::vector<double> const &points, size_t numPoints) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_ABSTRACTPOINTORDERING_HPP_ */
