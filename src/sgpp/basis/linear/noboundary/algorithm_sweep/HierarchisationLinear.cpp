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

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "basis/linear/noboundary/algorithm_sweep/HierarchisationLinear.hpp"

namespace sg
{

namespace detail
{

HierarchisationLinear::HierarchisationLinear(GridStorage* storage) : storage(storage)
{
}

HierarchisationLinear::~HierarchisationLinear()
{
}

void HierarchisationLinear::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	rec(source, result, index, dim, 0.0, 0.0);
}

void HierarchisationLinear::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{
	// current position on the grid
	size_t seq = index.seq();
	// value in the middle, needed for recursive call and calculation of the hierarchical surplus
	double fm = source[seq];

	// recursive calls for the right and left side of the current node
	if(index.hint() == false)
	{
		// descend left
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fl, fm);
		}

		// descend right
		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fm, fr);
		}

		// ascend
		index.up(dim);
	}

	// hierarchisation
	result[seq] = fm - ((fl + fr)/2.0);
}

}	// namespace detail

}	// namespace sg
