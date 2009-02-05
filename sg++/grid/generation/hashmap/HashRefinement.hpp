/*****************************************************************************/
/* This file is part of sg++, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sg++ is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sg++ is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with Foobar; if not, write to the Free Software                     */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef HASHREFINEMENT_HPP
#define HASHREFINEMENT_HPP

#include "grid/GridStorage.hpp"
#include "grid/generation/RefinementFunctor.hpp"

#include "exception/generation_exception.hpp"

#include <vector>
#include <cmath>

namespace sg
{

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

}

#endif /* HASHREFINEMENT_HPP */
