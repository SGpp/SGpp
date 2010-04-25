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

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"

namespace sg
{

namespace detail
{

PhiPhiDownBBLinearBoundary::PhiPhiDownBBLinearBoundary(GridStorage* storage) : PhiPhiDownBBLinear(storage)
{
}

PhiPhiDownBBLinearBoundary::~PhiPhiDownBBLinearBoundary()
{
}

void PhiPhiDownBBLinearBoundary::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	double q = this->boundingBox->getIntervalWidth(dim);
	double t = this->boundingBox->getIntervalOffset(dim);

	bool useBB = false;

	if (q != 1.0 || t != 0.0)
	{
		useBB = true;
	}

	// get boundary values
	double left_boundary;
	double right_boundary;
	size_t seq_left;
	size_t seq_right;

	/*
	 * Handle Level 0
	 */
	// This handles the diagonal only
	//////////////////////////////////////
	// left boundary
	index.left_levelzero(dim);
	seq_left = index.seq();
	left_boundary = source[seq_left];

	// right boundary
	index.right_levelzero(dim);
	seq_right = index.seq();
	right_boundary = source[seq_right];

	if (useBB)
	{
		// check boundary conditions
		if (this->boundingBox->hasDirichletBoundaryLeft(dim))
		{
			result[seq_left] = 0.0; //left_boundary
		}
		else
		{
			result[seq_left] = ((1.0/3.0)*left_boundary)*q;
		}
		if (this->boundingBox->hasDirichletBoundaryRight(dim))
		{
			result[seq_right] = 0.0; //right_boundary;
		}
		else
		{
			result[seq_right] = ((1.0/3.0)*right_boundary)*q;

			// down
			//////////////////////////////////////
			result[seq_right] += ((1.0/6.0)*left_boundary)*q;
		}

		// move to root
		if (!index.hint())
		{
			index.top(dim);

			if(!this->storage->end(index.seq()))
			{
				recBB(source, result, index, dim, left_boundary, right_boundary, q, t);
			}

			index.left_levelzero(dim);
		}
	}
	else
	{
		// check boundary conditions
		if (this->boundingBox->hasDirichletBoundaryLeft(dim))
		{
			result[seq_left] = 0.0; //left_boundary
		}
		else
		{
			result[seq_left] = (1.0/3.0)*left_boundary;
		}
		if (this->boundingBox->hasDirichletBoundaryRight(dim))
		{
			result[seq_right] = 0.0; //right_boundary;
		}
		else
		{
			result[seq_right] = (1.0/3.0)*right_boundary;

			// down
			//////////////////////////////////////
			result[seq_right] += (1.0/6.0)*left_boundary;
		}

		// move to root
		if (!index.hint())
		{
			index.top(dim);

			if(!this->storage->end(index.seq()))
			{
				rec(source, result, index, dim, left_boundary, right_boundary);
			}

			index.left_levelzero(dim);
		}
	}
}

} // namespace detail

} // namespace sg
