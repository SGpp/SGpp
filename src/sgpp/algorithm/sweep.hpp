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

#ifndef SWEEP_HPP
#define SWEEP_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

#include <vector>
#include <utility>
#include <iostream>

namespace sg
{

/**
 * Standard sweep operation
 * FUNC should be a class with overwritten operator(). For an example see laplace_up_functor in laplace.hpp.
 * It must be default constructable or copyable.
 * STORAGE must provide a grid_iterator supporting left_child, step_right, up, hint and seq.
 */
template<class FUNC>
class sweep
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	FUNC functor;
	GridStorage* storage;

public:
	/**
	 * Create a new sweep object with a default constructed functor
	 */
	sweep(GridStorage* storage) : functor(), storage(storage)
	{
	}

	/**
	 * Create a new sweep object with a copied functor
	 */
	sweep(FUNC& functor, GridStorage* storage) : functor(functor), storage(storage)
	{
	}

	~sweep()
	{
	}

	/**
	 * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep
	 */
	void sweep1D(DataVector& source, DataVector& result, size_t dim_sweep)
	{
		// generate a list of all dimension (-dim_sweep) from dimension recursion unrolling
		std::vector<size_t> dim_list;
		for(size_t i = 0; i < storage->dim(); i++)
		{
			if(i != dim_sweep)
			{
				dim_list.push_back(i);
			}
		}

		grid_iterator index(storage);

		sweep_rec(source, result, index, dim_list, storage->dim()-1, dim_sweep);

	}

protected:

	void sweep_rec(DataVector& source, DataVector& result, grid_iterator& index,
				std::vector<size_t>& dim_list, size_t dim_rem, size_t dim_sweep)
	{
		functor(source, result, index, dim_sweep);

		// dimension recursion unrolled
		for(size_t d = 0; d < dim_rem; d++)
		{
			size_t current_dim = dim_list[d];

			if(index.hint(current_dim))
			{
				continue;
			}

			index.left_child(current_dim);
			if(!storage->end(index.seq()))
			{
				sweep_rec(source, result, index, dim_list, d+1, dim_sweep);
			}

			index.step_right(current_dim);
			if(!storage->end(index.seq()))
			{
				sweep_rec(source, result, index, dim_list, d+1, dim_sweep);
			}

			index.up(current_dim);
		}
	}

};

}

#endif /* SWEEP_HPP */
