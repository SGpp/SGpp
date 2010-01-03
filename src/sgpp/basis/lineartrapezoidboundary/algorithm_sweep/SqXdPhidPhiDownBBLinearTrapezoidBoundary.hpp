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

#ifndef SQXDPHIDPHIDOWNBBLINEARTRAPEZOIDBOUNDARY_HPP
#define SQXDPHIDPHIDOWNBBLINEARTRAPEZOIDBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * down-operation in dimension dim. for use with sweep
 */
class SqXdPhidPhiDownBBLinearTrapezoidBoundary
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the GridStorage Object
	GridStorage* storage;
	/// Pointer to the bounding box Obejct
	BoundingBox* boundingBox;
	/// width of the interval in dimension
	double q;
	/// intervals offset in dimension
	double t;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	SqXdPhidPhiDownBBLinearTrapezoidBoundary(GridStorage* storage) : storage(storage), boundingBox(storage->getBoundingBox()), q(1.0), t(0.0)
	{
	}

	/**
	 * Destructor
	 */
	~SqXdPhidPhiDownBBLinearTrapezoidBoundary()
	{
	}

	/**
	 * This operations performs the calculation of down in the direction of dimension <i>dim</i>
	 *
	 * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
	 * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
	 * result)
	 * So please assure that both functions do exist!
	 *
	 * On level zero the getfixDirechletBoundaries of the storage object evaluated
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the down operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		q = boundingBox->getIntervalWidth(dim);
		t = boundingBox->getIntervalOffset(dim);

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
			double bbFactor = ((q*q) + (3.0*q*t) + (3.0*t*t))/(q);

			// check boundary conditions
			if (boundingBox->hasDirichletBoundaryLeft(dim))
			{
				result[seq_left] = 0.0; //left_boundary;
			}
			else
			{
				result[seq_left] = (1.0/3.0)*left_boundary*bbFactor;
			}

			if (boundingBox->hasDirichletBoundaryRight(dim))
			{
				result[seq_right] = 0.0; //right_boundary;
			}
			else
			{
				result[seq_right] = (1.0/3.0)*right_boundary*bbFactor;
				// down
				//////////////////////////////////////
				result[seq_right] -= (1.0/3.0)*left_boundary*bbFactor;
			}

			// move to root
			if (!index.hint())
			{
				index.top(dim);

				if(!storage->end(index.seq()))
				{
					recBB(source, result, index, dim, left_boundary, right_boundary);
				}

				index.left_levelzero(dim);
			}
		}
		else
		{
			// check boundary conditions
			if (boundingBox->hasDirichletBoundaryLeft(dim))
			{
				result[seq_left] = 0.0; //left_boundary;
			}
			else
			{
				result[seq_left] = (1.0/3.0)*left_boundary;
			}

			if (boundingBox->hasDirichletBoundaryRight(dim))
			{
				result[seq_right] = 0.0; //right_boundary;
			}
			else
			{
				result[seq_right] = (1.0/3.0)*right_boundary;
				// down
				//////////////////////////////////////
				result[seq_right] -= (1.0/3.0)*left_boundary;
			}

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
	}

protected:

	/**
	 * recursive function for the calculation of Down without Bounding Box support
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param fl function value on the left boundary
	 * @param fr function value on the right boundary
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
	{
		size_t seq = index.seq();

		double alpha_value = source[seq];

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double i_dbl = static_cast<double>(i);
		int l_int = static_cast<int>(l);

		double diagonal = ((1.0/3.0) + (i_dbl*i_dbl))*pow(2.0, 1-l_int);

		// integration
		result[seq] = (  (((1.0/pow(2.0, l_int))* i_dbl) * (fl-fr)) + (diagonal * alpha_value) );

		// dehierarchisation
		double fm = ((fl+fr)/2.0) + alpha_value;

		if(!index.hint())
		{
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fl, fm);
			}

			index.step_right(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fm, fr);
			}

			index.up(dim);
		}
	}

	/**
	 * recursive function for the calculation of Down with Bounding Box support
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param fl function value on the left boundary
	 * @param fr function value on the right boundary
	 */
	void recBB(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
	{
		size_t seq = index.seq();

		double alpha_value = source[seq];

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double i_dbl = static_cast<double>(i);
		int l_int = static_cast<int>(l);

		double diagonal = (1.0/3.0) * ((((pow(2.0, (1-l_int)))*q*q)*(3.0*(i_dbl*i_dbl) + 1)) + (12.0*t*q*i_dbl) + (3.0*(pow(2.0, (1+l_int)))*t*t))/(q);

		// integration
		result[seq] = (  (((1.0/pow(2.0, l_int))* i_dbl*q + t) * (fl-fr)) + (diagonal * alpha_value) );

		// dehierarchisation
		double fm = ((fl+fr)/2.0) + alpha_value;

		if(!index.hint())
		{
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				recBB(source, result, index, dim, fl, fm);
			}

			index.step_right(dim);
			if(!storage->end(index.seq()))
			{
				recBB(source, result, index, dim, fm, fr);
			}

			index.up(dim);
		}
	}
};

} // namespace detail

} // namespace sg

#endif /* SQXDPHIDPHIDOWNBBLINEARTRAPEZOIDBOUNDARY_HPP */
