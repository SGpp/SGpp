/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
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

#ifndef DEHIERARCHISATIONLINEARBOUNDARY_HPP
#define DEHIERARCHISATIONLINEARBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "basis/linear/noboundary/algorithm_sweep/DehierarchisationLinear.hpp"

namespace sg
{

namespace detail
{

/**
 * Class that implements the dehierarchisation on a linear sparse grid with boundaries. Therefore
 * the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal
 */
class DehierarchisationLinearBoundary : public DehierarchisationLinear
{
public:
	/**
	 * Constructor, must be bind to a grid
	 *
	 * @param storage the grid storage object of the the grid, on which the dehierarchisation should be executed
	 */
	DehierarchisationLinearBoundary(GridStorage* storage) : DehierarchisationLinear(storage)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~DehierarchisationLinearBoundary()
	{
	}

	/**
	 * Implements operator() needed by the sweep class during the grid traversal. This function
	 * is applied to the whole grid.
	 *
	 * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
	 * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
	 * result)
	 * So please assure that both functions do exist!
	 *
	 * @param source this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
	 * @param result this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		double left_boundary;
		double right_boundary;
		size_t seq;

		// left boundary
		index.left_levelzero(dim);
		seq = index.seq();
		left_boundary = source[seq];
		// right boundary
		index.right_levelzero(dim);
		seq = index.seq();
		right_boundary = source[seq];

		// move to root
		if (!index.hint())
		{
			index.top(dim);

			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, left_boundary, right_boundary);
			}

			index.left_levelzero(dim);
		}
	}
};

}	// namespace detail

}	// namespace sg

#endif /* DEHIERARCHISATIONLINEARBOUNDARY_HPP */
