/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007-2009 Dirk Pflueger (dirk.pflueger@in.tum.de)           */
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

#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

namespace sg
{

namespace detail
{

PhiPhiUpBBLinearBoundary::PhiPhiUpBBLinearBoundary(GridStorage* storage) : PhiPhiUpBBLinear(storage)
{
}

PhiPhiUpBBLinearBoundary::~PhiPhiUpBBLinearBoundary()
{
}

void PhiPhiUpBBLinearBoundary::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	double q = this->boundingBox->getIntervalWidth(dim);
	double t = this->boundingBox->getIntervalOffset(dim);

	bool useBB = false;

	if (q != 1.0 || t != 0.0)
	{
		useBB = true;
	}

	// get boundary values
	double fl = 0.0;
	double fr = 0.0;

	if (useBB)
	{
		if(!index.hint())
		{
			index.top(dim);

			if(!this->storage->end(index.seq()))
			{
				recBB(source, result, index, dim, fl, fr, q, t);
			}

			index.left_levelzero(dim);
		}

		size_t seq_left;
		size_t seq_right;

		// left boundary
		seq_left = index.seq();

		// right boundary
		index.right_levelzero(dim);
		seq_right = index.seq();

		// up
		//////////////////////////////////////
		//Left
		if (this->boundingBox->hasDirichletBoundaryLeft(dim))
		{
			result[seq_left] = 0.0; // source[seq_left];
		}
		else
		{
			result[seq_left] = fl;
			result[seq_left] += (((1.0/6.0)*source[seq_right])*q);
		}
		// Right
		if (this->boundingBox->hasDirichletBoundaryRight(dim))
		{
			result[seq_right] = 0.0; //source[seq_right];
		}
		else
		{
			result[seq_right] = fr;
		}

		index.left_levelzero(dim);
	}
	else
	{
		if(!index.hint())
		{
			index.top(dim);

			if(!this->storage->end(index.seq()))
			{
				rec(source, result, index, dim, fl, fr);
			}

			index.left_levelzero(dim);
		}

		size_t seq_left;
		size_t seq_right;

		// left boundary
		seq_left = index.seq();

		// right boundary
		index.right_levelzero(dim);
		seq_right = index.seq();

		// up
		//////////////////////////////////////
		//Left
		if (this->boundingBox->hasDirichletBoundaryLeft(dim))
		{
			result[seq_left] = 0.0; // source[seq_left];
		}
		else
		{
			result[seq_left] = fl;
			result[seq_left] += ((1.0/6.0)*source[seq_right]);
		}
		// Right
		if (this->boundingBox->hasDirichletBoundaryRight(dim))
		{
			result[seq_right] = 0.0; //source[seq_right];
		}
		else
		{
			result[seq_right] = fr;
		}

		index.left_levelzero(dim);
	}
}

} // namespace detail

} // namespace sg
