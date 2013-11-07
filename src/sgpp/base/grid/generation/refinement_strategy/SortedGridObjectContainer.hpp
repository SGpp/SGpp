/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)

//errorContainer object that stores predicted error for each Child of refinable points/subspaces.

#ifndef SORTEDGRIDOBJECTCONTAINER_HPP_
#define SORTEDGRIDOBJECTCONTAINER_HPP_

#include <stdlib.h>
#include <map>
#include <stdexcept>
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/refinement_strategy/RefinementDecorator.hpp"

using namespace std;

namespace sg {
namespace base {


class gridPointErrorContainer
{
public:
	AbstractRefinement::index_type* gridPoint;

	RefinementFunctor::value_type& operator [] (int i)
	{
		switch (i) {
			case 0: return leftChildErrorEstimate;
				break;
			case 1: return rightChildErrorEstimate;
				break;
			default:
					throw out_of_range("Invalid Index. This object only has allowed indices 0 and 1");
		}
	}

	RefinementFunctor::value_type largestErrorValue()
	{
		return max(leftChildErrorEstimate,rightChildErrorEstimate);
	}

	int largestErrorIndex()
	{
		return leftChildErrorEstimate < rightChildErrorEstimate ? 0 : 1;
	}

private:

	RefinementFunctor::value_type leftChildErrorEstimate;
	RefinementFunctor::value_type rightChildErrorEstimate;


};


// compares the errors of gridPointContainers. An errorContainer is larger then another, if the errorEstimate
// for one of his children is larger (e.g. if Poicont A has errors 1;3 and B has 0;2 A is larger).
// The compare operation returns true if the first argument is smaller then the second.
struct compareGridPointErrorContainer
{
public:

	//returns true if the first argument is smaller then the second.
	bool operator ()(const RefinementFunctor::value_type& keyA, const RefinementFunctor::value_type& keyB)
	{

		return keyA < keyB;
	}

};

typedef multimap<RefinementFunctor::value_type,gridPointErrorContainer,compareGridPointErrorContainer> GridObjectsSortedByError;
typedef pair<RefinementFunctor::value_type, gridPointErrorContainer> GridObjectErrorObject;


class SortedGridObjectContainer {
public:
	SortedGridObjectContainer();
	virtual ~SortedGridObjectContainer();

	void removeAddedErrorIndicators(size_t refinements_num,
									RefinementFunctor::value_type* addedMaxErrorValues,
									RefinementDecorator::index_type* addedMaxErrorIndices);

	void collectHighestErrorIndicators(size_t refinementes_num,
									   RefinementFunctor::value_type* maxErrorValues,
									   RefinementDecorator::index_type* maxErrorIndices);

private:
	GridObjectsSortedByError gridObjectsByError;

	/**
	 * Returns the index of the first accurance of minimal element in array.
	 * Used to find which entry is to be replaced next searching the maximum ones.
	 *
	 * @param array array with values
	 * @param length length of array
	 *
	 * @return index of the first occurrence of minimal element in array
	 */
	size_t getIndexOfMin(RefinementFunctor::value_type* array, size_t length) {
		size_t min_idx = 0;

		for (size_t i = 1; i < length; i++) {
			if (array[i] < array[min_idx])
				min_idx = i;
		}

		return min_idx;
	}

};

} /* namespace base */
} /* namespace sg */
#endif /* SORTEDGRIDOBJECTCONTAINER_HPP_ */
