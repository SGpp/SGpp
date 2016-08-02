/*
 * IterationPolicy.hpp
 *
 *  Created on: 12.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ITERATIONPOLICY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ITERATIONPOLICY_HPP_

#include "../common/AbstractPermutationIterator.hpp"
#include <vector>
#include <cstddef>
#include <memory>

namespace SGPP {
namespace combigrid {

/**
 * With the IterationPolicy class, an AbstractMultiStorage can be configured to traverse its entries in a different order than the identity order.
 * Therefore, an IterationPolicy can store a AbstractPermutationIterator for each dimension, which, if available, may define an other order in this dimension.
 */
class IterationPolicy {
	std::vector<std::shared_ptr<AbstractPermutationIterator>> iterators;

public:
	IterationPolicy() :
			iterators() {

	}

	/**
	 * @return Returns true iff there is a custom iterator in the given dimension.
	 */
	bool hasCustomIterator(size_t d) const {
		return d < iterators.size() && iterators[d];
	}

	/**
	 * @return Returns an iterator if there is one and nullptr if there is no custom iterator for this dimension.
	 */
	std::shared_ptr<AbstractPermutationIterator> getIterator(size_t d) const {
		if (d >= iterators.size()) {
			return nullptr;
		}

		return iterators[d];
	}

	/**
	 * Stores a custom iterator in the given dimension.
	 */
	void setIterator(size_t d, std::shared_ptr<AbstractPermutationIterator> iterator) {
		if (d >= iterators.size()) {
			iterators.resize(d + 1);
		}

		iterators[d] = iterator;
	}

	/**
	 * Advances the iterator in dimension d, if available, to the next position.
	 */
	void moveToNext(size_t d) {
		if (d < iterators.size()) {
			auto &it = iterators[d];
			if (it) {
				it->moveToNext();
			}
		}
	}

	/**
	 * Advances the iterator in dimension d, if available, to the next position and returns the permuted index at this position.
	 * If there is no custom iterator in dimension d, it returns the default value.
	 */
	size_t moveAndGetValue(size_t d, size_t defaultValue) {
		if (d < iterators.size()) {
			auto &it = iterators[d];
			if (it) {
				it->moveToNext();
				return it->value();
			}
		}

		return defaultValue;
	}

	/**
	 * @return Returns the value in dimension d if there is a custom iterator, or else the default value.
	 */
	size_t value(size_t d, size_t defaultValue) {
		if (d < iterators.size()) {
			auto &it = iterators[d];
			if (it) {
				return it->value();
			}
		}

		return defaultValue;
	}

	/**
	 * Resets the iterator in dimension d and returns its value, if available. Returns the default value otherwise.
	 */
	size_t resetAndGetValue(size_t d, size_t defaultValue) {
		if (d < iterators.size()) {
			auto &it = iterators[d];
			if (it) {
				it->reset();
				return it->value();
			}
		}

		return defaultValue;
	}

	/**
	 * Resets the iterator in dimension d to the start, if available.
	 */
	void reset(size_t d) {
		if (d < iterators.size()) {
			auto &it = iterators[d];
			if (it) {
				it->reset();
			}
		}
	}

	/**
	 * Default Iteration policy (no custom iterators).
	 */
	static IterationPolicy Default;
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ITERATIONPOLICY_HPP_ */
