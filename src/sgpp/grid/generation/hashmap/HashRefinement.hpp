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
/* along with sg++; if not, write to the Free Software                       */
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
 * Standard free refinement class for sparse grids with and without boundaries
 */
class HashRefinement
{
public:
	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;

	/// specifies, whether the refinement is executed on a grid with boundaries
	bool bWithBoundaries;

	/**
	 * Returns the index of the first accurance of minimal element in array
	 *
	 * @param: array array with ???
	 *
	 * @return: index of the first accurance of minimal element in array
	 *
	 * @todo: check if this has to be public
	 */
	int getIndexOfMin(RefinementFunctor::value_type* array)
	{
		int length = sizeof(array)/sizeof(RefinementFunctor::value_type);
		int min_idx = 0;
		for (int i = 1; i < length; i++)
		{
			if(array[i] < array[min_idx])
				min_idx = i;
		}

		return min_idx;
	}

	/**
	 * Performs the refinement on grid
	 *
	 * @param: storage hashmap that stores the grid points
	 * @param: functor a function used to determine if refinement is needed
	 * @param: bBoundaries specifies, whether the grid has boundaries or not
	 */
	void free_refine(GridStorage* storage, RefinementFunctor* functor, bool bBoundaries = false)
	{
		if(storage->size() == 0)
		{
			throw generation_exception("storage empty");
		}

		//Algorithm should be able to look for several points in grid to refine
		//So we store an array with refinements_num maximal points
		int refinements_num = functor->getRefinementsNum();
		RefinementFunctor::value_type max_values[refinements_num];
		size_t max_indexes[refinements_num];
		for (int i = 0; i<refinements_num; i++){
			max_values[i] = functor->start();
			max_indexes[i] = 0;
		}
		int min_idx = 0;

		RefinementFunctor::value_type max_value = max_values[min_idx];
		size_t max_index = max_indexes[min_idx];

		bWithBoundaries = bBoundaries;

		index_type index;
		GridStorage::grid_map_iterator end_iter = storage->end();

		// do a lot of copy&paste due to performance reasons
		if (bWithBoundaries == true)
		{
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

					if (source_level == 0)
					{
						// we only have one child on level 1
						index.set(d, 1, 1);
						child_iter = storage->find(&index);
						// if there no more grid points --> test if we should refine the grid
						if(child_iter == end_iter)
						{
							RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
							if(current_value > max_value)
							{
								//Replace the minimal point in result array, find the new  minimal point
								max_values[min_idx] = current_value;
								max_indexes[min_idx] = iter->second;
								min_idx = getIndexOfMin(max_values);
								max_value = max_values[min_idx];
								break;
							}
						}
					}
					else
					{
						// left child
						index.set(d, source_level + 1, 2 * source_index - 1);
						child_iter = storage->find(&index);
						// if there no more grid points --> test if we should refine the grid
						if(child_iter == end_iter)
						{
							RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
							if(current_value > max_value)
							{
								//Replace the minimal point in result array, find the new  minimal point
								max_values[min_idx] = current_value;
								max_indexes[min_idx] = iter->second;
								min_idx = getIndexOfMin(max_values);
								max_value = max_values[min_idx];
								break;
							}
						}

						// right child
						index.set(d, source_level + 1, 2 * source_index + 1);
						child_iter = storage->find(&index);
						if(child_iter == end_iter)
						{
							RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
							if(current_value > max_value)
							{
								//Replace the minimal point in result array, find the new minimal point
								max_values[min_idx] = current_value;
								max_indexes[min_idx] = iter->second;
								min_idx = getIndexOfMin(max_values);
								max_value = max_values[min_idx];
								break;
							}
						}
					}

					index.set(d, source_level, source_index);
				}
			}
		}
		else
		{
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

					// left child
					index.set(d, source_level + 1, 2 * source_index - 1);
					child_iter = storage->find(&index);
					// if there no more grid points --> test if we should refine the grid
					if(child_iter == end_iter)
					{
						RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
						if(current_value > max_value)
						{
							//Replace the minimal point in result array, find the new  minimal point
							max_values[min_idx] = current_value;
							max_indexes[min_idx] = iter->second;
							min_idx = getIndexOfMin(max_values);
							max_value = max_values[min_idx];
							break;
						}
					}

					// right child
					index.set(d, source_level + 1, 2 * source_index + 1);
					child_iter = storage->find(&index);
					if(child_iter == end_iter)
					{
						RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
						if(current_value > max_value)
						{
							//Replace the minimal point in result array, find the new minimal point
							max_values[min_idx] = current_value;
							max_indexes[min_idx] = iter->second;
							min_idx = getIndexOfMin(max_values);
							max_value = max_values[min_idx];
							break;
						}
					}

					index.set(d, source_level, source_index);
				}
			}
		}

		//can refine grid on several points
		for (int i = 0; i < refinements_num; i++){
			max_value = max_values[i];
			max_index = max_indexes[i];
			if(max_value > functor->start())
			{
				refine_gridpoint(storage, max_index);
			}
		}

	}


protected:
	/**
	 * This method refines a grid point be generating the children in every dimension
	 * of the grid.
	 *
	 * @param: storage hashmap that stores the gridpoints
	 * @param: refine_index the index in the hashmap of the point that should be refined
	 */
	void refine_gridpoint(GridStorage* storage, size_t refine_index)
	{
		index_type index((*storage)[refine_index]);

		// @todo: Maybe it's possible to move predecessor/successor discovery into the storage concept
		for(size_t d = 0; d < storage->dim(); d++)
		{
			index_t source_index;
			level_t source_level;
			index.get(d, source_level, source_index);

			if (source_level == 0)
			{
				// we only have one child on level 1
				index.set(d, 1, 1);
				if(!storage->has_key(&index))
				{
					create_gridpoint(storage, index);
				}
			}
			else
			{
				// generate left child, if necessary
				index.set(d, source_level + 1, 2 * source_index - 1);
				if(!storage->has_key(&index))
				{
					create_gridpoint(storage, index);
				}

				// generate right child, if necessary
				index.set(d, source_level + 1, 2 * source_index + 1);
				if(!storage->has_key(&index))
				{
					create_gridpoint(storage, index);
				}
			}

			index.set(d, source_level, source_index);
		}
	}

	/**
	 * This method creates a new point on the grid. It checks if some parents or
	 * children are needed in other dimensions.
	 *
	 * @param storage hashmap that stores the gridpoinrs
	 * @param index the point that should be inserted
	 *
	 * @todo 	check, if it's better the implement an own class for grid refinement on a grid
	 * 			with boundaries.
	 */
	void create_gridpoint(GridStorage* storage, index_type& index)
	{
		for(size_t d = 0; d < storage->dim(); d++)
		{
			index_t source_index;
			level_t source_level;
			index.get(d, source_level, source_index);

			if (source_level == 1)
			{
				// check if we need some additional points on the boundaries, only needed on a N dim grid
				if (bWithBoundaries == true && storage->dim() > 1)
				{
					// test if there are boundaries in every dimension for this grid point
					// left boundary
					index.set(d, 0, 0);
					if(!storage->has_key(&index))
					{
						create_gridpoint(storage, index);
					}

					// right boundary
					index.set(d, 0, 1);
					if(!storage->has_key(&index))
					{
						create_gridpoint(storage, index);
					}

					// restore values
					index.set(d, source_level, source_index);
				}
			}

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

				// restore values
				index.set(d, source_level, source_index);
			}
		}
		storage->insert(index);
	}

};

}

#endif /* HASHREFINEMENT_HPP */
