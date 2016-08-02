/*
 * SortedPermutationIterator.hpp
 *
 *  Created on: 18.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_SORTEDPERMUTATIONITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_SORTEDPERMUTATIONITERATOR_HPP_

#include <sgpp/combigrid/common/AbstractPermutationIterator.hpp>
#include <vector>

namespace SGPP {
namespace combigrid {

class SortedPermutationIterator: public AbstractPermutationIterator {
	std::vector<size_t> permutation;
	size_t currentIndex;

	SortedPermutationIterator(std::vector<size_t> const &permutation, size_t currentIndex);
public:
	SortedPermutationIterator(std::vector<SGPP::float_t> const &points, size_t numPoints);
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
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_SORTEDPERMUTATIONITERATOR_HPP_ */
