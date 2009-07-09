/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef HIERARCHISATIONMODLINEAR_HPP
#define HIERARCHISATIONMODLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * Class that implements the hierarchisation on a mod linear sparse grid. Therefore
 * the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal
 */
class HierarchisationModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// the grid object
	GridStorage* storage;

public:
	/**
	 * Constructor, must be bind to a grid
	 *
	 * @param storage the grid storage object of the the grid, on which the hierarchisation should be executed
	 */
	HierarchisationModLinear(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~HierarchisationModLinear()
	{
	}

	/**
	 * Implements operator() needed by the sweep class during the grid traversal. This function
	 * is applied to the whole grid.
	 *
	 * @param source this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
	 * @param result this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		rec(source, result, index, dim, 0.0, 0.0);
	}

protected:

	/**
	 * Recursive hierarchisaton algorithm, this algorithms works in-place -> source should be equal to result
	 *
	 * @todo (heinecke, nice) add graphical explanation here
	 *
	 * @param source this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
	 * @param result this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 * @param fl left value of the current region regarded in this step of the recursion
	 * @param fr right value of the current region regarded in this step of the recursion
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
	{
		// current position on the grid
		size_t seq = index.seq();
		// value in the middle, needed for recursive call and calculation of the hierarchical surplus
		double fm = source[seq];

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		// recursive calls for the right and left side of the current node
		if(index.hint() == false)
		{
			double fltemp = fl;
			double frtemp = fr;

			// When we descend the hierarchical basis we have to modify the boundary values
			// in case the index is 1 or (2^l)-1 or we are on the first level
			// level 1, constant function
			if(l == 1)
			{
				// constant function
				fltemp = fm;
				frtemp = fm;
			}
			// left boundary
			else if(i == 1)
			{
				double ftemp;
				ftemp = fr - fm;
				fltemp = fm - ftemp;
			}
			// right boundary
			else if(static_cast<int>(i) == static_cast<int>((1 << l)-1))
			{
				double ftemp;
				ftemp = fl - fm;
				frtemp = fm - ftemp;
			}
			// inner functions
			else
			{
			}

			// descend left
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fltemp, fm);
			}

			// descend right
			index.step_right(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fm, frtemp);
			}

			// ascend
			index.up(dim);
		}

		// hierarchisation
		result[seq] = fm - ((fl + fr)/2.0);
	}
};

}	// namespace detail

}	// namespace sg

#endif /* HIERARCHISATIONMODLINEAR_HPP */
