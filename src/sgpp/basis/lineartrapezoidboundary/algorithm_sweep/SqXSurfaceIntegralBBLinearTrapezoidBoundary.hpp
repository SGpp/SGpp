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

#ifndef SQXSURFACEINTEGRALBBLINEARTRAPEZOIDBOUNDARY_HPP
#define SQXSURFACEINTEGRALBBLINEARTRAPEZOIDBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * x^2 Surface Integration in dimension dim
 */
class SqXSurfaceIntegralBBLinearTrapezoidBoundary
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Pointer to the bounding box Obejct
	BoundingBox* boundingBox;
	/// width of the interval in dimension
	double q;
	/// intervals offset in dimension = left boundary value
	double t;
	/// right boundary value = t + q
	double r;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	SqXSurfaceIntegralBBLinearTrapezoidBoundary(GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()), q(1.0), t(0.0), r(1.0)
	{
	}

	/**
	 * Destructor
	 */
	~SqXSurfaceIntegralBBLinearTrapezoidBoundary()
	{
	}

	/**
	 * This operations performs the calculation the surface integral in dimension <i>dim</i>
	 *
	 * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
	 * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
	 * result)
	 * So please assure that both functions do exist!
	 *
	 * On level zero the getfixDirechletBoundaries of the storage object evaluated
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the up operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		q = boundingBox->getIntervalWidth(dim);
		t = boundingBox->getIntervalOffset(dim);
		r = t+q;

		size_t seq_left;
		size_t seq_right;

		/*
		 * Handle Level 0
		 */
		// left boundary
		index.left_levelzero(dim);
		seq_left = index.seq();

		// right boundary
		index.right_levelzero(dim);
		seq_right = index.seq();

		//check boundary conditions
		if (boundingBox->hasDirichletBoundaryLeft(dim))
		{
			result[seq_left] = 0.0; //source[seq_left];
		}
		else
		{
			result[seq_left] = (-1.0)*(source[seq_left]*(t*t));
		}

		if (boundingBox->hasDirichletBoundaryRight(dim))
		{
			result[seq_right] = 0.0; //source[seq_right];
		}
		else
		{
			result[seq_right] = (source[seq_right]*(r*r));
		}

		// move to root
		if (!index.hint())
		{
			index.top(dim);

			if(!storage->end(index.seq()))
			{
				// set everything else to zero
				rec(source, result, index, dim);
			}

			index.left_levelzero(dim);
		}
	}

protected:

	/**
	 * recursive function for the calculation of the surface integral
	 *
	 * traverse all basis function in dimension <i>dim</i> with level greater than zero
	 * and set their coefficients to zero.
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		size_t seq = index.seq();

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		result[seq] = 0.0;

		if(!index.hint())
		{
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim);
			}

			index.step_right(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim);
			}

			index.up(dim);
		}

	}
};

} // namespace detail

} // namespace sg

#endif /* SQXSURFACEINTEGRALBBLINEARTRAPEZOIDBOUNDARY_HPP */
