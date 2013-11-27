/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "PredictiveRefinement.hpp"

namespace sg {
namespace base {


void PredictiveRefinement::free_refineTest(GridStorage* storage,
		PredictiveRefinementIndicator* errorIndicator)
{
	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	// the functor->getRefinementsNum() largest grid points should be refined.
	// gather them in an array max_values
	size_t refinements_num = errorIndicator->getRefinementsNum();
	// values
	RefinementFunctor::value_type* max_values = new RefinementFunctor::value_type[refinements_num];
	// indices
	size_t* max_indices = new size_t [refinements_num];

	// initialization
	for (size_t i = 0; i < refinements_num; i++) {
		max_values[i] = errorIndicator->start();
		max_indices[i] = 0;
	}

	collectRefinablePoints(storage, errorIndicator, refinements_num, max_indices, max_values);
	// now refine all grid points which satisfy the refinement criteria
	refineGridpointsCollection(storage, errorIndicator, refinements_num, max_indices, max_values);
	delete [] max_values;
	delete [] max_indices;
}

void PredictiveRefinement::collectRefinablePoints(
		GridStorage* storage, PredictiveRefinementIndicator* errorIndicator,
		size_t refinements_num, size_t* max_indices,
		RefinementFunctor::value_type* max_values)
{
	size_t min_idx = 0;

		// max value equals min value
		RefinementFunctor::value_type max_value = max_values[min_idx];
		//size_t max_index = max_indices[min_idx];

		index_type index;
		GridStorage::grid_map_iterator end_iter = storage->end();

		// start iterating over whole grid
		for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
			index = *(iter->first);

			GridStorage::grid_map_iterator child_iter;

			// @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
			// check for each grid point whether it can be refined (i.e., whether not all kids exist yet)
			// if yes, check whether it belongs to the refinements_num largest ones
			for (size_t d = 0; d < storage->dim(); d++) {
				index_t source_index;
				level_t source_level;
				index.get(d, source_level, source_index);

				std::cout << "analyzing point " << index.toString() << "\n";

				// test existence of left child
				index.set(d, source_level + 1, 2 * source_index - 1);
				std::cout << "testing existence of left child " << index.toString() << "\n";
				child_iter = storage->find(&index);

				// if there no more grid points --> test if we should refine the grid
				if (child_iter == end_iter) {
					RefinementFunctor::value_type current_value = (*errorIndicator)(&index);

					std::cout << "is refineable\n";

					if (current_value > max_value) {
						// replace the minimal point in result array, find the new  minimal point
						max_values[min_idx] = current_value;
						max_indices[min_idx] = iter->second;
						min_idx = getIndexOfMin(max_values, refinements_num);
						max_value = max_values[min_idx];
						break;
					}
				}

				// test existance of right child
				index.set(d, source_level + 1, 2 * source_index + 1);
				std::cout << "testing existence of right child " << index.toString() << "\n";
				child_iter = storage->find(&index);

				if (child_iter == end_iter) {
					RefinementFunctor::value_type current_value = (*errorIndicator)(&index);
					std::cout << "is refineable\n";

					if (current_value > max_value) {
						// replace the minimal point in result array, find the new minimal point
						max_values[min_idx] = current_value;
						max_indices[min_idx] = iter->second;
						min_idx = getIndexOfMin(max_values, refinements_num);
						max_value = max_values[min_idx];
						break;
					}
				}

				// reset current grid point in dimension d
				index.set(d, source_level, source_index);
			}
		}
}

void PredictiveRefinement::refineGridpointsCollection(
		GridStorage* storage,
		PredictiveRefinementIndicator* errorIndicator,
		size_t refinements_num, size_t* max_indices,
		RefinementFunctor::value_type* max_values)
{
	RefinementFunctor::value_type max_value;
	size_t max_index;
	// now refine all grid points which satisfy the refinement criteria
	double threshold = errorIndicator->getRefinementThreshold();

	for (size_t i = 0; i < refinements_num; i++) {
		max_value = max_values[i];
		max_index = max_indices[i];

		if (max_value > errorIndicator->start() && fabs(max_value) >= threshold) {
			refineGridpoint(storage, max_index);
		}
	}
}

} /* namespace base */
} /* namespace sg */
