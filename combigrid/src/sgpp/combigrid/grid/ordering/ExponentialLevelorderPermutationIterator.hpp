/*
 * ExponentialLevelorderPermutationIterator.hpp
 *
 *  Created on: 19.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPERMUTATIONITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPERMUTATIONITERATOR_HPP_

#include <sgpp/combigrid/common/AbstractPermutationIterator.hpp>

namespace sgpp{
namespace combigrid {

class ExponentialLevelorderPermutationIterator: public AbstractPermutationIterator {
	size_t currentIndex;
	size_t level;
	size_t numPoints;

public:
	ExponentialLevelorderPermutationIterator(size_t level, size_t numPoints, size_t currentIndex = 0);

	virtual ~ExponentialLevelorderPermutationIterator();

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

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_GRID_POINTS_ORDERING_EXPONENTIALLEVELORDERPERMUTATIONITERATOR_HPP_ */
