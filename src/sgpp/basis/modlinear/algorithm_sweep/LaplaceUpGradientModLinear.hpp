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

#ifndef LAPLACEUPGRADIENTMODLINEAR_HPP
#define LAPLACEUPGRADIENTMODLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

namespace sg
{

namespace detail
{

/**
 * Implements the upGradient Method needed for the Laplace operator on mod linear grids
 */
class LaplaceUpGradientModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	/// Pointer to GridStorage object
	GridStorage* storage;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	LaplaceUpGradientModLinear(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~LaplaceUpGradientModLinear()
	{
	}

	/**
	 * This operations performs the calculation of upGradient in the direction of dimension <i>dim</i>
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the down operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		double f = 0.0;
		rec(source, result, index, dim, f);
	}

protected:

	/**
	 * recursive function for the calculation of upGradient
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param f function value in the middle
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& f)
	{
		size_t seq = index.seq();

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double alpha_value = source[seq];
		double ht = pow(2.0, l);

		if(l == 1)
		{
			f = 0.0;
			if(!index.hint())
			{
				index.left_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}

				f = 0.0;
				index.step_right(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}

			result[seq] = 0.0;
		}
		// left boundary
		else if(i == 1)
		{
			f = 0.0;
			if(!index.hint())
			{
				index.left_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}

			result[seq] = ht * f;

			f += 2.0 * alpha_value;
		}
		// right boundary
		else if(static_cast<int>(i) == static_cast<int>((1 << l)-1))
		{
			f = 0.0;
			if(!index.hint())
			{
				index.right_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}

			result[seq] = ht * f;

			f += 2.0 * alpha_value;
		}
	}

};

}

}

#endif /* LAPLACEUPGRADIENTMODLINEAR_HPP */
