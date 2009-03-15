/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

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
 * @todo add picture here
 *
 * 2. a modified boundary grid with level 0 and a pentagon cut
 * trough the sub space scheme.
 *
 * @todo add picture here
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
	 * @todo level should be of type level_t but swig doesnt want that
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
			index.push(d, 1, 1);
		}

		this->regular_rec(storage, index, storage->dim() - 1, storage->dim(), level + storage->dim() - 1);
	}

	/**
	 * Generates a regular sparse grid of level levels with boundaries
	 *
	 * @todo: level should be of type level_t but swig doesnt want that
	 *
	 * @param storage Hashmap, that stores the grid points
	 * @param level maximum level of the sparse grid
	 * @param UScaledBoundaries true -> generate sparse grid with less points on the boundary, pentagon cut through subspace scheme
	 */
	void regularWithBoundaries(GridStorage* storage, int level, bool UScaledBoundaries)
	{
		if(storage->size() > 0)
		{
			throw generation_exception("storage not empty");
		}

		index_type index(storage->dim());

		if (UScaledBoundaries == true)
		{
			if (level == 0)
			{
				for(size_t d = 0; d < storage->dim(); d++)
				{
					index.push(d, 0, 0);
				}

				this->boundaries_rec(storage, index, storage->dim() - 1, 0, 0);
			}
			else
			{
				for(size_t d = 0; d < storage->dim(); d++)
				{
					index.push(d, 1, 1);
				}

				this->boundaries_UScaled_rec(storage, index, storage->dim()-1, storage->dim(), level + storage->dim() - 1);
			}
		}
		else
		{
			/* new grid generation
			 *
			 * for all level the same calculation of the level sum is implemented:
			 * |l| = n + d -1
			 */
			for(size_t d = 0; d < storage->dim(); d++)
			{
				index.push(d, 0, 0);
			}

			this->boundaries_rec(storage, index, storage->dim()-1, 0, level);

			/*
			 * old grid generation
			 *
			 * here a switching in calculating the level sum
			 * is implemented: For level 0 |l| = n and for all other levels |l| = n+d-1
			 * are applied.
			 */
			/*if (level < 2)
			{
				for(size_t d = 0; d < storage->dim(); d++)
				{
					index.push(d, 0, 0);
				}

				this->boundaries_rec(storage, index, storage->dim()-1, 0, level);
			}
			else
			{
				for(size_t d = 0; d < storage->dim(); d++)
				{
					index.push(d, 1, 1);
				}

				// handle 1 D grid
				if (storage->dim() == 1)
				{
					this->boundariesFull_rec(storage, index, storage->dim()-1, storage->dim(), level + storage->dim() - 1, true, false);
				}
				else
				{
					this->boundariesFull_rec(storage, index, storage->dim()-1, storage->dim(), level + storage->dim() - 1, false, false);
				}
			}*/
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
            for(index_t i = 1; i <= 1<<(l-1); i++)
            {
                index.push(0, l, 2*i-1);
                storage->insert(index);
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
	 */
	void boundaries_UScaled_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level)
	{
		if(current_dim == 0)
		{
			boundaries_UScaled_rec_1d(storage, index, current_level, level);
		}
		else
		{
			index_t source_index;
			level_t source_level;

			index.get(current_dim, source_level, source_index);

			if(current_level <= level)
			{
				if (source_level == 1)
				{
					index.push(current_dim, 0, 0);
					this->boundaries_UScaled_rec(storage, index, current_dim-1, current_level, level);

					index.push(current_dim, 0, 1);
					this->boundaries_UScaled_rec(storage, index, current_dim-1, current_level, level);

					index.push(current_dim, source_level, source_index);
				}

				// d-1 recursion
				this->boundaries_UScaled_rec(storage, index, current_dim - 1, current_level, level);
			}

			if(current_level < level)
			{
				index.push(current_dim, source_level + 1, 2*source_index - 1);
				this->boundaries_UScaled_rec(storage, index, current_dim, current_level + 1, level);

				index.push(current_dim, source_level + 1, 2*source_index + 1);
				this->boundaries_UScaled_rec(storage, index, current_dim, current_level + 1, level);
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
	 */
	void boundaries_UScaled_rec_1d(GridStorage* storage, index_type& index, level_t current_level, level_t level)
	{
		for(level_t l = 0; l <= level-current_level + 1; l++)
		{
			if (l == 0)
			{
				index.push(0, 0, 0);
				storage->insert(index);
				index.push(0, 0, 1);
				storage->insert(index);
			}
			else
			{
				for(index_t i = 1; i <= 1<<(l-1); i++)
				{
					index.push(0, l, 2*i-1);
					storage->insert(index);
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
			// d-1 recursion
			if (source_level == 0)
			{
				if (current_dim == 0)
				{
					index.push(0, 0, 0);
					storage->insert(index);
					index.push(0, 0, 1);
					storage->insert(index);

					index.push(current_dim, source_level, source_index);
				}
				else
				{
					this->boundaries_rec(storage, index, current_dim-1, current_level, level);

					index.push(current_dim, 0, 1);
					this->boundaries_rec(storage, index, current_dim-1, current_level, level);

					index.push(current_dim, source_level, source_index);
				}
			}
			else
			{
				if (current_dim == 0)
				{
					storage->insert(index);
				}
				else
				{
					this->boundaries_rec(storage, index, current_dim-1, current_level, level);
				}
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

	/*
	 * recursive construction of the spare grid with boundaries, classic level 0 approach, only for level greater 1
	 *
	 * @param storage hashmap that stores the grid points
	 * @param index point's index
	 * @param current_dim current working dimension
	 * @param current_level current level in this construction step
	 * @param level maximum level of the sparse grid
	 */
	/*
	void boundariesFull_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level, bool tatooine, bool kessel)
	{
		if(current_dim == 0)
		{
			boundariesFull_rec_lastd(storage, index, current_level, level, tatooine, kessel);
		}
		else
		{
			index_t source_index;
			level_t source_level;

			index.get(current_dim, source_level, source_index);

			if(current_level <= level)
			{
				if (source_level == 1)
				{
					index.push(current_dim, 0, 0);
					this->boundariesFull_rec(storage, index, current_dim-1, current_level, level, true, kessel);

					index.push(current_dim, 0, 1);
					this->boundariesFull_rec(storage, index, current_dim-1, current_level, level, true, kessel);

					index.push(current_dim, source_level, source_index);
				}

				// d-1 recursion
				this->boundariesFull_rec(storage, index, current_dim - 1, current_level, level, tatooine, kessel);
			}

			if(current_level < level)
			{
				if (current_level == (level-1))
				{
					index.push(current_dim, source_level + 1, 2*source_index - 1);
					this->boundariesFull_rec(storage, index, current_dim, current_level + 1, level, tatooine, true);

					index.push(current_dim, source_level + 1, 2*source_index + 1);
					this->boundariesFull_rec(storage, index, current_dim, current_level + 1, level, tatooine, true);
				}
				else
				{
					index.push(current_dim, source_level + 1, 2*source_index - 1);
					this->boundariesFull_rec(storage, index, current_dim, current_level + 1, level, tatooine, kessel);

					index.push(current_dim, source_level + 1, 2*source_index + 1);
					this->boundariesFull_rec(storage, index, current_dim, current_level + 1, level, tatooine, kessel);
				}
			}

			index.push(current_dim, source_level, source_index);
		}
	}*/

	/*
	 * generate points of the last dimension (dim == 0), version of the classic sparse grid
	 * with level 0, overall level > 1
	 *
	 * @param storage the hashmap that stores the grid points
	 * @param index point's index that should be created on the grid
	 * @param current_level current level of the grid generation
	 * @param level maximum level of grid
	 */
	/*
	void boundariesFull_rec_lastd(GridStorage* storage, index_type& index, level_t current_level, level_t level, bool tatooine, bool kessel)
	{
		if (tatooine == false && kessel == true)
		{
			index.push(0, 0, 0);
			storage->insert(index);
			index.push(0, 0, 1);
			storage->insert(index);
		}
		else
		{
			if (tatooine == false)
				level--;

			for(level_t l = 0; l <= level-current_level + 1; l++)
			{
				if (l == 0)
				{
					index.push(0, 0, 0);
					storage->insert(index);
					index.push(0, 0, 1);
					storage->insert(index);
				}
				else
				{
					for(index_t i = 1; i <= 1<<(l-1); i++)
					{
						index.push(0, l, 2*i-1);
						storage->insert(index);
					}
				}
			}
		}
	}*/
};

}

#endif /* HASHGENERATOR_HPP */
