/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2007  Jörg Blank (blankj@in.tum.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GENERATION_HPP_
#define GENERATION_HPP_

#include "DataVector.h"
#include "GridStorage.hpp"
#include "Operations.hpp"
#include "exceptions.hpp"

#include <vector>
#include <cmath>

namespace sg
{


/**
 * This abstracts the refinement criteria out of the refinement algorithm
 */
class SurplusRefinementFunctor : public RefinementFunctor
{
public:


	SurplusRefinementFunctor(DataVector* alpha) : alpha(alpha)
	{
	}
	
	virtual ~SurplusRefinementFunctor() {}
	
	/**
	 * This should be returning a refinement value for every grid point.
	 * The point with the highest value will be refined.
	 */
	virtual double operator()(GridStorage* storage, size_t seq)
	{
		return fabs(alpha->get(seq));
	}
	
	/**
	 * This should return a bottom.
	 */
	virtual double start()
	{
		return 0.0;
	}
	
	
protected:
	DataVector* alpha;
};


/**
 * Standard free refinement class
 */
class HashRefinement
{
public:
	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;
	
	
	void free_refine(GridStorage* storage, RefinementFunctor* functor)
	{
		if(storage->size() == 0)
		{
			throw generation_exception("storage empty");
		}
		
		RefinementFunctor::value_type max_value = functor->start();
		size_t max_index = 0;
		
		index_type index;
		GridStorage::grid_map_iterator end_iter = storage->end();
		
		// I think this may be depedent on local support
		for(GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++)
		{
			index = *(iter->first);
			
			GridStorage::grid_map_iterator child_iter;
			
			// TODO: Maybe it's possible to move predecessor/successor discovery into the storage concept
			for(size_t d = 0; d < storage->dim(); d++)
			{
				index_t source_index;
				level_t source_level;
				index.get(d, source_level, source_index);
				
				index.set(d, source_level + 1, 2 * source_index - 1);
				child_iter = storage->find(&index);			
				if(child_iter == end_iter)
				{
					RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
					if(current_value > max_value)
					{
						max_value = current_value;
						max_index = iter->second;
						break;
					} 
				}
				
				index.set(d, source_level + 1, 2 * source_index + 1);
				child_iter = storage->find(&index);			
				if(child_iter == end_iter)
				{
					RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
					if(current_value > max_value)
					{
						max_value = current_value;
						max_index = iter->second;
						break;
					} 
				}
				
				index.set(d, source_level, source_index);
			}
		}
		
		if(max_value > functor->start())
		{
			refine_gridpoint(storage, max_index);
		}		
		
	}
	
	
protected:
	
	void refine_gridpoint(GridStorage* storage, size_t refine_index)
	{
		index_type index((*storage)[refine_index]);
		
		// TODO: Maybe it's possible to move predecessor/successor discovery into the storage concept
		for(size_t d = 0; d < storage->dim(); d++)
		{
			index_t source_index;
			level_t source_level;
			index.get(d, source_level, source_index);
			
			index.set(d, source_level + 1, 2 * source_index - 1);
			if(!storage->has_key(&index))
			{
				create_gridpoint(storage, index);
			}

			index.set(d, source_level + 1, 2 * source_index + 1);
			if(!storage->has_key(&index))
			{
				create_gridpoint(storage, index);
			}

			index.set(d, source_level, source_index);
		}
	}
	
	void create_gridpoint(GridStorage* storage, index_type& index)
	{
		for(size_t d = 0; d < storage->dim(); d++)
		{
			index_t source_index;
			level_t source_level;
			index.get(d, source_level, source_index);
		
			if(source_level > 1)
			{
				// TODO: Maybe it's possible to move predecessor/successor discovery into the storage concept
				if(((source_index + 1) / 2) % 2 == 1)
				{
					index.set(d, source_level - 1, (source_index + 1) / 2);
				}
				else
				{
					index.set(d, source_level - 1, (source_index - 1) / 2);
				}
				
				if(!storage->has_key(&index))
				{
					create_gridpoint(storage, index);	
				}
				index.set(d, source_level, source_index);
			}
		}		
		storage->insert(index);
	}
	
};


class HashGenerator
{
public:

	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;


/*
 * Generates a regular sparse grid of level levels
 * TODO: level should be of type level_t but swig doesnt want that
 */
	void regular(GridStorage* storage, int level)
	{
		if(storage->size() > 0)
		{
			throw generation_exception("storage not empty");
		}
		
		index_type index(storage->dim());
		
		for(size_t d = 0; d < storage->dim(); d++)
		{
			index.push(d, 1, 1);
		}
		
		this->regular_rec(storage, index, storage->dim() - 1, storage->dim(), level + storage->dim() - 1);
	}
	
	
	
protected:
	void regular_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level) 
	{
		if(current_dim == 0)
		{
			regular_rec_1d(storage, index, current_level, level);
		}
		else
		{
			index_t source_index;
			level_t source_level;
			
			index.get(current_dim, source_level, source_index);
			
			if(current_level <= level)
			{
				// d-1 recursion
				this->regular_rec(storage, index, current_dim - 1, current_level, level);
			}
			
			if(current_level < level)
			{
				// current_level + 1 recursion
				index.push(current_dim, source_level + 1, 2*source_index - 1); 
				this->regular_rec(storage, index, current_dim, current_level + 1, level);
				
				index.push(current_dim, source_level + 1, 2*source_index + 1);
				this->regular_rec(storage, index, current_dim, current_level + 1, level);
				
			}
			
			index.push(current_dim, source_level, source_index);
		}
	}
	
	void regular_rec_1d(GridStorage* storage, index_type& index, level_t current_level, level_t level)
	{
        for(level_t l = 1; l <= level-current_level + 1; l++)
        {
            for(index_t i = 1; i <= 1<<(l-1); i++)
            {
                index.push(0, l, 2*i-1);
                storage->insert(index);
            }
        }		
	}
};
	
}


#endif /*GENERATION_HPP_*/
