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
/* GNU General Public License for more details.                              */
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

class HashGenerator
{
public:

	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;


	/**
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

	/**
	 * Generates a regular sparse grid of level levels with boundaries
	 * TODO: level should be of type level_t but swig doesnt want that
	 */
	void regularWithBoundaries(GridStorage* storage, int level, bool fullBoundaries)
	{
		if(storage->size() > 0)
		{
			throw generation_exception("storage not empty");
		}

		index_type index(storage->dim());

		// the boundaries are hyperplanes
		for(size_t d = 0; d < storage->dim(); d++)
		{
			index.push(d, 1, 1);
		}

		if (fullBoundaries == true)
		{
			throw new generation_exception("full boundaries are not supported, yet!");
			//this->boundariesFull_rec(storage, index, storage->dim() - 1, storage->dim(), level + storage->dim() - 1, level + storage->dim());
		}
		else
		{
			this->boundaries_rec(storage, index, storage->dim() - 1, storage->dim(), level + storage->dim() - 1);
		}
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

	void boundaries_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level)
	{
		if(current_dim == 0)
		{
			boundaries_rec_1d(storage, index, current_level, level);
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
					this->boundaries_rec(storage, index, current_dim-1, current_level, level);

					index.push(current_dim, 0, 1);
					this->boundaries_rec(storage, index, current_dim-1, current_level, level);

					index.push(current_dim, source_level, source_index);
				}

				// d-1 recursion
				this->boundaries_rec(storage, index, current_dim - 1, current_level, level);
			}

			if(current_level < level)
			{
				index.push(current_dim, source_level + 1, 2*source_index - 1);
				this->boundaries_rec(storage, index, current_dim, current_level + 1, level);

				index.push(current_dim, source_level + 1, 2*source_index + 1);
				this->boundaries_rec(storage, index, current_dim, current_level + 1, level);
			}

			index.push(current_dim, source_level, source_index);
		}
	}

	void boundariesFull_rec(GridStorage* storage, index_type& index, size_t current_dim, level_t current_level, level_t level, level_t level_boundary)
	{
		if(current_dim == 0)
		{
			boundaries_rec_1d(storage, index, current_level, level);
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
					this->boundariesFull_rec(storage, index, current_dim-1, current_level, level_boundary, level_boundary);

					index.push(current_dim, 0, 1);
					this->boundariesFull_rec(storage, index, current_dim-1, current_level, level_boundary, level_boundary);

					index.push(current_dim, source_level, source_index);
				}

				// d-1 recursion
				this->boundariesFull_rec(storage, index, current_dim - 1, current_level, level, level_boundary);
			}

			if(current_level < level)
			{
				index.push(current_dim, source_level + 1, 2*source_index - 1);
				this->boundariesFull_rec(storage, index, current_dim, current_level + 1, level, level_boundary);

				index.push(current_dim, source_level + 1, 2*source_index + 1);
				this->boundariesFull_rec(storage, index, current_dim, current_level + 1, level, level_boundary);
			}

			index.push(current_dim, source_level, source_index);
		}
	}


	void boundaries_rec_1d(GridStorage* storage, index_type& index, level_t current_level, level_t level)
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

	void boundariesFull_rec_lastd(GridStorage* storage, index_type& index, level_t current_level, level_t level)
	{
		for(level_t l = 1; l <= level-current_level + 1; l++)
		{
			for(index_t i = 1; i <= 1<<(l-1); i++)
			{
				index.push(storage->dim()-1, l, 2*i-1);
				storage->insert(index);
			}
		}
	}
};

}

#endif /* HASHGENERATOR_HPP */
