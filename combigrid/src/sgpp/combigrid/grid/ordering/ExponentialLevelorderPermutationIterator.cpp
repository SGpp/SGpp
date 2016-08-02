/*
 * ExponentialLevelorderPermutationIterator.cpp
 *
 *  Created on: 19.12.2015
 *      Author: david
 */

#include "ExponentialLevelorderPermutationIterator.hpp"

namespace SGPP {
namespace combigrid {

ExponentialLevelorderPermutationIterator::ExponentialLevelorderPermutationIterator(size_t level, size_t numPoints, size_t currentIndex)
	: currentIndex(currentIndex)
	, level(level)
	, numPoints(numPoints) {
}

ExponentialLevelorderPermutationIterator::~ExponentialLevelorderPermutationIterator() {
}

void ExponentialLevelorderPermutationIterator::reset() {
	currentIndex = 0;
}

size_t ExponentialLevelorderPermutationIterator::value() {
	size_t lastPoint = numPoints - 1;

	if(level == 0) {
		return 0;
	}

	if(currentIndex == 0) {
		return 1;
	}
	if(currentIndex == lastPoint) {
		return 2;
	}
	if(currentIndex*2 == lastPoint) {
		return 0; // this is the midpoint of the interval
	}

	// We're at a level >= 2, so we can pretend that the boundary points are the first two points and the rest is just like without boundary
	size_t currentLevel = level;
	size_t dividedIndex = currentIndex;

	// loop is not endless because we checked for currentIndex == 0 already
	while(dividedIndex % 2 == 0) {
		dividedIndex /= 2;
		--currentLevel;
	}

	// the divided index is always odd, so we can "compress" it
	dividedIndex = (dividedIndex - 1)/2;

	// two boundary points, offset 2^l - 1 (= sum_{k=0}^{l-1} 2^k), dividedIndex starts from 1, so we have to subtract 1
	return 2 + ((1L << (currentLevel-1)) - 1) + dividedIndex;
}

void ExponentialLevelorderPermutationIterator::moveToNext() {
	++currentIndex;
}

std::shared_ptr<AbstractPermutationIterator> ExponentialLevelorderPermutationIterator::clone() {
	return std::shared_ptr<AbstractPermutationIterator>(new ExponentialLevelorderPermutationIterator(level, numPoints, currentIndex));
}

} /* namespace combigrid */
} /* namespace SGPP */
