/*
 * ExponentialLevelorderPointOrdering.hpp
 *
 *  Created on: 19.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPOINTORDERING_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPOINTORDERING_HPP_

#include "AbstractPointOrdering.hpp"

namespace sgpp{
namespace combigrid {

class ExponentialLevelorderPointOrdering: public AbstractPointOrdering {
public:
	virtual ~ExponentialLevelorderPointOrdering();

	virtual size_t convertIndex(size_t level, size_t numPoints, size_t index);

	virtual size_t numPoints(size_t level);

	virtual std::shared_ptr<AbstractPermutationIterator> getSortedPermutationIterator(size_t level, std::vector<double> const &points,
			size_t numPoints);
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPOINTORDERING_HPP_ */
