/*
 * IdentityPointOrdering.hpp
 *
 *  Created on: 18.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_IDENTITYPOINTORDERING_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_IDENTITYPOINTORDERING_HPP_

#include "AbstractPointOrdering.hpp"
#include "../growth/AbstractGrowthStrategy.hpp"


namespace SGPP {
namespace combigrid {

class IdentityPointOrdering: public AbstractPointOrdering {
	std::shared_ptr<AbstractGrowthStrategy> growthStrategy;
	bool isSorted;

public:
	IdentityPointOrdering(std::shared_ptr<AbstractGrowthStrategy> growthStrategy, bool isSorted);

	virtual ~IdentityPointOrdering();

	virtual size_t convertIndex(size_t level, size_t numPoints, size_t index);

	virtual size_t numPoints(size_t level);

	virtual std::shared_ptr<AbstractPermutationIterator> getSortedPermutationIterator(size_t level, std::vector<SGPP::float_t> const &points,
			size_t numPoints);
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_IDENTITYPOINTORDERING_HPP_ */
