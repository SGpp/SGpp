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

#include "basis/linear/boundary/algorithm_sweep/DehierarchisationLinearBoundary.hpp"

namespace sg
{

namespace detail
{

DehierarchisationLinearBoundary::DehierarchisationLinearBoundary(GridStorage* storage) : DehierarchisationLinear(storage)
{
}

DehierarchisationLinearBoundary::~DehierarchisationLinearBoundary()
{
}

void DehierarchisationLinearBoundary::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
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

}	// namespace detail

}	// namespace sg
