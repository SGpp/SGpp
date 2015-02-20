// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include "GSGRefinement.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

 void GSGRefinement::collectRefinablePoints(
		GridStorage* storage, RefinementFunctor* functor,
		size_t refinements_num, size_t* max_indices,
		RefinementFunctor::value_type* max_values)
{
	//this refinement algorithm uses the predictive refinement indicator.
	//dynamic casting is used to maintain the signature of the algorithm, but still be able to use the
	//predictive refinement indicator with it.
	//unused
//	 PredictiveRefinementIndicator* errorIndicator = dynamic_cast<PredictiveRefinementIndicator*>(functor);

	 //unused
//	size_t min_idx = 0;

	// max value equals min value
	//unused
//	RefinementFunctor::value_type max_value = max_values[min_idx];

	index_type index;
	//unused
//	GridStorage::grid_map_iterator end_iter = storage->end();

	// TODO Valeriy: find less intrusive mechanism to handle the error contribution
//	// start iterating over whole grid
//	for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
//		index = *(iter->first);
//
//
//
//		GridStorage::grid_map_iterator child_iter;
//		ErrorType newPointError(index);
//
//		// @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
//		// check for each grid point whether it can be refined (i.e., whether not all kids exist yet)
//		// if yes, check whether it belongs to the refinements_num largest ones
//		for (size_t d = 0; d < storage->dim(); d++) {
//			index_t source_index;
//			level_t source_level;
//			index.get(d, source_level, source_index);
//
//			// test existence of left child
//			index.set(d, source_level + 1, 2 * source_index - 1);
//			child_iter = storage->find(&index);
//
//			if (child_iter == end_iter) {
//				//use the predictive error indicator, which takes a pointer to the grid point object
//				//instead of the storage index
//				newPointError += (*errorIndicator)(&index);
//			}
//
//			// test existance of right child
//			index.set(d, source_level + 1, 2 * source_index + 1);
//			child_iter = storage->find(&index);
//
//			if (child_iter == end_iter) {//use predictive refinement indicator
//				//use the predictive error indicator, which takes a pointer to the grid point object
//				//instead of the storage index
//				newPointError += (*errorIndicator)(&index);
//			}
//			// reset current grid point in dimension d
//			index.set(d, source_level, source_index);
//		}
//
//		if (newPointError > max_value)
//		{
//			// replace the minimal point in result array, find the new minimal point
//			max_values[min_idx] = newPointError.getContribPerPoint();
//			max_indices[min_idx] = iter->second;
//			min_idx = getIndexOfMin(max_values, refinements_num);
//			max_value = max_values[min_idx];
//
//		}
//	}

	//DEBUG
//	std::cout << "selected points\n";
//	for(size_t i = 0; i < refinements_num;++i)
//	{
//		std::cout << "point: " << storage->get(max_indices[i])->toString()  << "with error:" << max_values[i] << "\n";
//	}

}

} /* namespace base */
} /* namespace SGPP */
