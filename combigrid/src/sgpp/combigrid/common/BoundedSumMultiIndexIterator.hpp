/*
 * BoundedSumMultiIndexIterator.hpp
 *
 *  Created on: 14.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_BOUNDEDSUMMULTIINDEXITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_BOUNDEDSUMMULTIINDEXITERATOR_HPP_

#include "../definitions.hpp"
#include <cstddef>

namespace sgpp{
namespace combigrid {

class BoundedSumMultiIndexIterator {
	MultiIndex index;
	size_t maxIndexSum;
	size_t indexSum;
	bool valid;

public:
	BoundedSumMultiIndexIterator(size_t dim, size_t maxIndexSum) : index(dim, 0), maxIndexSum(maxIndexSum), indexSum(0), valid(true) {

	}

	void reset() {
		valid = true;
		for(size_t i = 0; i < index.size(); ++i) {
			index[i] = 0;
		}
	}

	bool isValid() {
		return valid;
	}

	MultiIndex &value() {
		return index;
	}

	size_t indexAt(size_t d) const {
		return index[d];
	}

	int moveToNext() {
		size_t lastDim = index.size() - 1;
		size_t d = lastDim;

		while (indexSum >= maxIndexSum) {
			if (d == 0) {
				valid = false;
				return -1;
			}

			indexSum -= index[d];
			index[d] = 0;

			--d;
		}

		++index[d];
		++indexSum;

		return static_cast<int>(lastDim - d);
	}
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_BOUNDEDSUMMULTIINDEXITERATOR_HPP_ */
