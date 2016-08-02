/*
 * AbstractPointHierarchy.hpp
 *
 *  Created on: 04.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTHIERARCHY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTHIERARCHY_HPP_

#include <sgpp/globaldef.hpp>
#include <cstddef>
#include <sgpp/combigrid/common/AbstractPermutationIterator.hpp>

#include <vector>

namespace sgpp{
namespace combigrid {

/**
 * Contains all the necessary information about a family of one-dimensional grids:
 * - Point distribution
 * - Level-numPoints-mapping
 * - Point sorting
 * - Hierarchy information
 *
 * Derived classes of AbstractPointHierarchy also save the grid points, avoiding expensive re-computation.
 */
class AbstractPointHierarchy {
public:
	virtual ~AbstractPointHierarchy();

	/**
	 * @return Returns the grid point for the given level and index. (0 <= index < getNumPoints(level))
	 */
	virtual double getPoint(size_t level, size_t index) = 0;

	/**
	 * @return Returns a vector with all points for the given level. If sorted == true, then the points returned are sorted.
	 * Depending on the configured PointOrdering, this might be slower than not sorting the points.
	 */
	virtual std::vector<double> getPoints(size_t level, bool sorted) = 0;

	/**
	 * @return Returns the number of points in the given level.
	 */
	virtual size_t getNumPoints(size_t level) = 0;

	/**
	 * @return Returns true if the points of a level are always a subset of the points at the next level.
	 * This should also mean that the points that are common to two levels share the same indices.
	 */
	virtual bool isNested() = 0;

	/**
	 * @return Returns a permutation iterator which can be used to traverse the points in their sorted (ascending) order.
	 */
	virtual std::shared_ptr<AbstractPermutationIterator> getSortedPermutationIterator(size_t level) = 0;

};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ABSTRACTPOINTHIERARCHY_HPP_ */
