/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHGENERATOR_HPP
#define HASHGENERATOR_HPP

#include "grid/GridStorage.hpp"

#include "exception/generation_exception.hpp"

#include <vector>
#include <cmath>

namespace sg
{

/**
 * This class provides the generation functionality of sparse grids
 * based on hashmaps.
 *
 * Grids with and without boundaries are supported. For boundary grids
 * two cases are supported:
 *
 * 1. the classic sparse grid with level 0 and a diagonal
 * cut through the sub space scheme.
 *
 * @todo (heinecke, nice) add picture here
 *
 * 2. a modified boundary grid with level 0 and a pentagon cut
 * trough the sub space scheme.
 *
 * @todo (heinecke, nice) add picture here
 */
class HashGenerator
{
public:

	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;


	/**
	 * Generates a regular sparse grid of level levels, without boundaries
	 *
	 * @todo (blank) level should be of type level_t but swig doesnt want that
	 *
	 * @param storage Hashmap, that stores the grid points
	 * @param level maximum level of the grid
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
			index.push(d, 1, 1, false);
		}

		size_t mydim = storage->dim();
		this->regular_rec(storage, index, ((storage->dim()) - 1), mydim, ((level + (storage->dim())) - 1));
	}

	/**
	 * Generates a regular sparse grid of level levels with boundaries
	 *
	 * @todo (blank) level should be of type level_t but swig doesnt want that
	 *
	 * @param storage Hashmap, that stores the grid points
	 * @param level maximum level of the sparse grid
	 * @param bTrapezoidBoundaries true -> generate sparse grid with less points on the boundary, pentagon cut through subspace scheme
	 */
	void regularWithBoundaries(GridStorage* storage, int level, bool bTrapezoidBoundaries)
	{
		if(storage->size() > 0)
		{
			throw generation_exception("storage not empty");
		}

		index_type index(storage->dim());

		if (bTrapezoidBoundaries == true)
		{
			if (level == 0)
			{
				for(size_t d = 0; d < storage->dim(); d++)
				{
					index.push(d, 0, 0, false);
				}

				this->boundaries_rec(storage, index, storage->dim() - 1, 0, 0);
			}
			else
			{
				for(size_t d = 0; d < storage->dim(); d++)
				{
					index.push(d, 1, 1, false);
				}

				this->boundaries_trapezoid_rec(storage, index, storage->dim()-1, storage->dim(), level + storage->dim() - 1, false);
			}
		}
		else
		{
			/* new grid generation
			 *
			 * for all level the same calculation of the level sum is implemented:
			 * |l| <= n
			 */
			for(size_t d = 0; d < storage->dim(); d++)
			{
				index.push(d, 0, 0, false);
			}

			this->boundaries_rec(storage, index, storage->dim()-1, 0, level);
		}
	}

protected:
	/**
	 * recursive construction of the spare grid without boundaries
	 *
	 * @param storage hashmap that stores the grid points
	 * @param index point's index
	 * @param current_dim current working dimension
	 * @param current_level current level in this construction step
	 * @param level maximum level of the sparse grid
	 */
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
				// set Leaf option of index
				if (current_level == level)
				{
					index.setLeaf(true);
				}
				else
				{
					index.setLeaf(false);
				}

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

	/**
	 * generate points of the last dimension (dim == 0), without boundaries
	 *
	 * @param storage the hashmap that stores the grid points
	 * @param index point's index that should be created on the grid
	 * @param current_level current level of the grid generation
	 * @param level maximum level of grid
	 */
	void regular_rec_1d(GridStorage* storage, index_type& index, level_t current_level, level_t level)
	{
        for(level_t l = 1; l <= level-current_level + 1; l++)
        {
            if (l == level-current_level+1)
            {
            	for(index_t i = 1; static_cast<int>(i) <= static_cast<int>(1<<(l-1)); i++)
                {
                    index.push(0, l, 2*i-1, true);
                    storage->insert(index);
                }
            }
            else
            {
            	for(index_t i = 1; static_cast<int>(i) <= static_cast<int>(1<<(l-1)); i++)
                {
                    index.push(0, l, 2*i-1, false);
                    storage->insert(index);
                }
            }
        }
	}

	/**
	 * recursive construction of the spare grid with boundaries, pentagon cut
	 *
	 * @param storage hashmap that stores the grid points
	 * @param index point's index
	 * @param current_dim current working dimension
	 * @param current_level current level in this construction step
	 * @param level maximum level of the sparse grid
	 * @param bLevelZero specifies if the current index has a level zero component
	 */
	void boundaries_trapezoid_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level, bool bLevelZero)
	{
		if(current_dim == 0)
		{
			boundaries_Trapezoid_rec_1d(storage, index, current_level, level, bLevelZero);
		}
		else
		{
			index_t source_index;
			level_t source_level;

			index.get(current_dim, source_level, source_index);

			if(current_level <= level)
			{
				// set Leaf option of index
				bool bLeafProperty = false;
				if (current_level == level)
				{
					bLeafProperty = true;
				}
				else
				{
					bLeafProperty = false;
				}

				if (source_level == 1)
				{
					index.push(current_dim, 0, 0, false);
					this->boundaries_trapezoid_rec(storage, index, current_dim-1, current_level, level, true);

					index.push(current_dim, 0, 1, false);
					this->boundaries_trapezoid_rec(storage, index, current_dim-1, current_level, level, true);

					index.push(current_dim, source_level, source_index);
				}

				// d-1 recursion
				index.setLeaf(bLeafProperty);
				this->boundaries_trapezoid_rec(storage, index, current_dim - 1, current_level, level, bLevelZero);
			}

			if(current_level < level)
			{
				index.push(current_dim, source_level + 1, 2*source_index - 1);
				this->boundaries_trapezoid_rec(storage, index, current_dim, current_level + 1, level, bLevelZero);

				index.push(current_dim, source_level + 1, 2*source_index + 1);
				this->boundaries_trapezoid_rec(storage, index, current_dim, current_level + 1, level, bLevelZero);
			}

			index.push(current_dim, source_level, source_index);
		}
	}

	/**
	 * generate points of the last dimension (dim == 0), version of pentagon cut in
	 * sub space scheme
	 *
	 * @param storage the hashmap that stores the grid points
	 * @param index point's index that should be created on the grid
	 * @param current_level current level of the grid generation
	 * @param level maximum level of grid
	 * @param bLevelZero specifies if the current index has a level zero component
	 */
	void boundaries_Trapezoid_rec_1d(GridStorage* storage, index_type& index, level_t current_level, level_t level, bool bLevelZero)
	{
		bool bLevelGreaterZero = !bLevelZero;
		for(level_t l = 0; l <= level-current_level + 1; l++)
		{
			if (l == level-current_level+1)
			{
				if (l == 0)
				{
					index.push(0, 0, 0, false);
					storage->insert(index);
					index.push(0, 0, 1, false);
					storage->insert(index);
				}
				else
				{
					for(index_t i = 1; static_cast<int>(i) <= static_cast<int>(1<<(l-1)); i++)
					{
						index.push(0, l, 2*i-1, (true && bLevelGreaterZero));
						storage->insert(index);
					}
				}
			}
			else
			{
				if (l == 0)
				{
					index.push(0, 0, 0, false);
					storage->insert(index);
					index.push(0, 0, 1, false);
					storage->insert(index);
				}
				else
				{
					for(index_t i = 1; static_cast<int>(i) <= static_cast<int>(1<<(l-1)); i++)
					{
						index.push(0, l, 2*i-1, false);
						storage->insert(index);
					}
				}
			}
		}
	}

	/**
	 * recursive construction of the spare grid with boundaries, classic level 0 approach, only for level 0 and 1
	 *
	 * @param storage hashmap that stores the grid points
	 * @param index point's index
	 * @param current_dim current working dimension
	 * @param current_level current level in this construction step
	 * @param level maximum level of the sparse grid
	 */
	void boundaries_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level)
	{
		index_t source_index;
		level_t source_level;

		index.get(current_dim, source_level, source_index);

		if(current_level <= level)
		{
			// set Leaf option of index
			bool bSaveLeafProperty = index.isLeaf();
			bool bLeafProperty = false;
			if (current_level == level)
			{
				bLeafProperty = true;
			}
			else
			{
				bLeafProperty = false;
			}

			// d-1 recursion
			if (source_level == 0)
			{
				if (current_dim == 0)
				{
					index.push(0, 0, 0, bLeafProperty);
					storage->insert(index);
					index.push(0, 0, 1, bLeafProperty);
					storage->insert(index);

					index.push(current_dim, source_level, source_index, bSaveLeafProperty);
				}
				else
				{
					index.push(current_dim, 0, 0, bLeafProperty);
					this->boundaries_rec(storage, index, current_dim-1, current_level, level);

					index.push(current_dim, 0, 1, bLeafProperty);
					this->boundaries_rec(storage, index, current_dim-1, current_level, level);

					index.push(current_dim, source_level, source_index, bSaveLeafProperty);
				}
			}
			else
			{
				index.setLeaf(bLeafProperty);
				if (current_dim == 0)
				{
					storage->insert(index);
				}
				else
				{
					this->boundaries_rec(storage, index, current_dim-1, current_level, level);
				}
				index.setLeaf(bSaveLeafProperty);
			}
		}

		if(current_level < level)
		{
			if (source_level == 0 && source_index == 0)
			{
				index.push(current_dim, source_level + 1, 1);
				this->boundaries_rec(storage, index, current_dim, current_level + 1, level);
			}
			else
			{
				index.push(current_dim, source_level + 1, 2*source_index - 1);
				this->boundaries_rec(storage, index, current_dim, current_level + 1, level);

				index.push(current_dim, source_level + 1, 2*source_index + 1);
				this->boundaries_rec(storage, index, current_dim, current_level + 1, level);
			}
		}

		index.push(current_dim, source_level, source_index);
	}
};

}

#endif /* HASHGENERATOR_HPP */
