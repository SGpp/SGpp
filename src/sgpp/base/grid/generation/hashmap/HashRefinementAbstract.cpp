/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy KHakhutskyy (khakhutv@in.tum.de)


#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"
#include "base/grid/generation/refinement_strategy/RefinementStrategy.hpp"


namespace sg
{
namespace base
{
size_t HashRefinementAbstract::getIndexOfMin(RefinementFunctor::value_type* array, size_t length)
{
		size_t min_idx = 0;
		for (size_t i = 1; i < length; i++)
		{
			if(array[i] < array[min_idx])
				min_idx = i;
		}

		return min_idx;
}

void HashRefinementAbstract::create_gridpoint_1d(index_type& index,
		size_t d, GridStorage * storage, index_t& source_index, level_t& source_level)
{
	index.get(d, source_level, source_index);
	if(source_level > 1){
		if(((source_index + 1) / 2) % 2 == 1){
			index.set(d, source_level - 1, (source_index + 1) / 2);
		}else{
			index.set(d, source_level - 1, (source_index - 1) / 2);
		}
		create_gridpoint_subroutine(storage, index);
		// restore values
		index.set(d, source_level, source_index);
	}
}

void HashRefinementAbstract::strategy_refine(GridStorage* storage,
										RefinementStrategy& refinement_strategy)
{
		refinement_strategy.refine(storage, this);
}

}
}


