/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "SortedGridObjectContainer.hpp"


namespace sg {
namespace base {

SortedGridObjectContainer::SortedGridObjectContainer() {


}

SortedGridObjectContainer::~SortedGridObjectContainer() {

}

void SortedGridObjectContainer::removeAddedErrorIndicators(
		size_t refinements_num,
		RefinementFunctor::value_type* addedMaxErrorValues,
		RefinementDecorator::index_type* addedMaxErrorIndices) {

	//iterate over the value array
	for (size_t i = 0; i < refinements_num; ++i) {

		//error is not zero => has been added. nullify error in
		if (addedMaxErrorValues[i]>0) {


			//find all objects with that specific error.
			pair<GridObjectsSortedByError::iterator,GridObjectsSortedByError::iterator> equalErrorObjects;
			equalErrorObjects = gridObjectsByError.equal_range(addedMaxErrorValues[i]);

			for (GridObjectsSortedByError::iterator it = equalErrorObjects.first;
					it!=equalErrorObjects.second;
					++it)
			{
				gridPointErrorContainer* storedValue = &(it->second);
				//check if it is the one we are after.

				//TODO umbau zu pointern!
				if(&(addedMaxErrorIndices[i]) == storedValue->gridPoint ){
					(*storedValue)[storedValue->largestErrorIndex()] = 0;

					if(storedValue->largestErrorValue()== 0){
						//delete element from multimap, if both indicators are zero.
						gridObjectsByError.erase(it);
					}else {
						//reformat key
						GridObjectErrorObject container;
						container.first=it->second.largestErrorValue();
						container.second = it->second;
						gridObjectsByError.erase(it);
						gridObjectsByError.insert(container);
					}
				}
			}
		}

	}
}

void SortedGridObjectContainer::collectHighestErrorIndicators(
		size_t refinementes_num, RefinementFunctor::value_type* maxErrorValues,
		RefinementDecorator::index_type* maxErrorIndices) {

	//reset error values;
	size_t smallestErrorIndex = 0;
	for (size_t i = 0; i < refinementes_num; ++i) {
		maxErrorValues = 0;
	}

	size_t numberOfAnalyizedGridPoints = 0;
	for(GridObjectsSortedByError::iterator it = gridObjectsByError.begin(); it!=gridObjectsByError.end(); ++it){
		if(numberOfAnalyizedGridPoints < refinementes_num){

			if (it->second[0]>maxErrorValues[smallestErrorIndex]) {
				maxErrorValues[smallestErrorIndex] = it->second[0];
				maxErrorIndices[smallestErrorIndex] = it->second.gridPoint;
				smallestErrorIndex = getIndexOfMin(maxErrorValues,refinementes_num);
			}
			if (it->second[1]>maxErrorValues[smallestErrorIndex]) {
				maxErrorValues[smallestErrorIndex] = it->second[1];
				maxErrorIndices[smallestErrorIndex] = it->second.gridPoint;
				smallestErrorIndex = getIndexOfMin(maxErrorValues,refinementes_num);
			}
			numberOfAnalyizedGridPoints++;
		}else {
			break;
		}
	}
}

} /* namespace base */
} /* namespace sg */
