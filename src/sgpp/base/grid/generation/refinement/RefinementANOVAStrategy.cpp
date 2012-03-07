/*
 * RefinementANOVAStrategy.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: khakhutv_local
 */

#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"
#include "base/grid/generation/refinement/RefinementANOVAStrategy.hpp"

namespace sg {
namespace base {

RefinementANOVAStrategy::RefinementANOVAStrategy() {
	// TODO Auto-generated constructor stub

}

RefinementANOVAStrategy::~RefinementANOVAStrategy() {
	// TODO Auto-generated destructor stub
}

void RefinementANOVAStrategy::refine(GridStorage* storage, HashRefinementAbstract* hash_refinement) {
	// only the function with the local error indicator grater than the
	// threshold will be refined
	double threshold = refinement_functor_->getRefinementThreshold();

	if(storage->size() == 0)
	{
		throw generation_exception("storage empty");
	}

	HashRefinementAbstract::index_type index;
	GridStorage::grid_map_iterator end_iter = storage->end();

	// start iterating over whole grid
	for(GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++)
	{
		index = *(iter->first);

		GridStorage::grid_map_iterator child_iter;

		for(size_t d = 0; d < storage->dim(); d++)
		{
			HashRefinementAbstract::index_t source_index;
			HashRefinementAbstract::level_t source_level;
			index.get(d, source_level, source_index);

			// in order to remain in the same ANOVA component, we shouldn't
			// refine the constant functions (level 1)
			if (source_level <= 1) continue;

			// test existence of the left child
			index.set(d, source_level + 1, 2 * source_index - 1);
			child_iter = storage->find(&index);
			// if there no more grid points --> test if we should refine the grid
			if(child_iter == end_iter)
			{
				RefinementFunctor::value_type current_value = fabs((*refinement_functor_)(storage, iter->second));
				if(current_value > threshold)
					refine_gridpoint_directional(storage, iter->second, d, hash_refinement);
				break; // you don't have to test for the right child
			}

			// if there is a left child test the existence of the right child
			index.set(d, source_level + 1, 2 * source_index + 1);
			child_iter = storage->find(&index);
			if(child_iter == end_iter)
			{
				RefinementFunctor::value_type current_value = fabs((*refinement_functor_)(storage, iter->second));
				if(current_value > threshold)
					refine_gridpoint_directional(storage, iter->second, d, hash_refinement);
				break;
			}

			// reset current grid point in dimension d
			index.set(d, source_level, source_index);
		}
	}
}

void RefinementANOVAStrategy::refine_gridpoint_directional(GridStorage* storage, size_t refine_index, size_t d, HashRefinementAbstract* hash_refinement)
	{
		HashRefinementAbstract::index_type index((*storage)[refine_index]);

		//Sets leaf property of index, which is refined to false
		(*storage)[refine_index]->setLeaf(false);


		HashRefinementAbstract::index_t source_index;
		HashRefinementAbstract::level_t source_level;
		index.get(d, source_level, source_index);

		// generate left child, if necessary
		index.set(d, source_level + 1, 2 * source_index - 1);
		if(!storage->has_key(&index))
		{
			index.setLeaf(true);
			hash_refinement->create_gridpoint(storage, index);
		}

		// generate right child, if necessary
		index.set(d, source_level + 1, 2 * source_index + 1);
		if(!storage->has_key(&index))
		{
			index.setLeaf(true);
			hash_refinement->create_gridpoint(storage, index);
		}

		index.set(d, source_level, source_index);

	}


} /* namespace base */
} /* namespace sg */
