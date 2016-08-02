/*
 * IdentityPointOrdering.cpp
 *
 *  Created on: 18.12.2015
 *      Author: david
 */

#include "IdentityPointOrdering.hpp"
#include "SortedPermutationIterator.hpp"

namespace SGPP {
namespace combigrid {

IdentityPointOrdering::IdentityPointOrdering(std::shared_ptr<AbstractGrowthStrategy> growthStrategy, bool isSorted)
	: growthStrategy(growthStrategy)
	, isSorted(isSorted) {
}

IdentityPointOrdering::~IdentityPointOrdering() {
}

size_t IdentityPointOrdering::convertIndex(size_t level, size_t numPoints, size_t index) {
	return index;
}

size_t IdentityPointOrdering::numPoints(size_t level) {
	return growthStrategy->numPoints(level);
}

std::shared_ptr<AbstractPermutationIterator> IdentityPointOrdering::getSortedPermutationIterator(size_t level,
		const std::vector<SGPP::float_t>& points, size_t numPoints) {
	if(isSorted) {
		return std::shared_ptr<AbstractPermutationIterator>(nullptr);
	}

	return std::shared_ptr<AbstractPermutationIterator>(new SortedPermutationIterator(points, numPoints));
}

} /* namespace combigrid */
} /* namespace SGPP */
