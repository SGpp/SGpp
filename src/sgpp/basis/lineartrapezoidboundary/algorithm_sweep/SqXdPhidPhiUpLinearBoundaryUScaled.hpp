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

#ifndef SQXDPHIDPHIUPLINEARBOUNDARYUSCALED_HPP
#define SQXDPHIDPHIUPLINEARBOUNDARYUSCALED_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * up-operation in dimension dim. for use with sweep
 */
class SqXdPhidPhiUpLinearBoundaryUScaled
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
	SqXdPhidPhiUpLinearBoundaryUScaled(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~SqXdPhidPhiUpLinearBoundaryUScaled()
	{
	}

	/**
	 * This operations performs the calculation of up in the direction of dimension <i>dim</i>
	 *
	 * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
	 * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
	 * result)
	 * So please assure that both functions do exist!
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the up operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		// get boundary values
		double fl = 0.0;
		double fr = 0.0;

		rec(source, result, index, dim, fl, fr);
	}

protected:

	/**
	 * recursive function for the calculation of Up
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 * @param fl function value on the left boundary, reference parameter
	 * @param fr function value on the right boundary, reference parameter
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr)
	{
		size_t seq = index.seq();

		fl = fr = 0.0;
		double fml = 0.0;
		double fmr = 0.0;

		GridStorage::index_type::level_type current_level;
		GridStorage::index_type::index_type current_index;

		index.get(dim, current_level, current_index);

		if(current_level > 0)
		{
			if(!index.hint())
			{
				index.left_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, fl, fml);
				}

				index.step_right(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, fmr, fr);
				}

				index.up(dim);
			}
		}
		else
		{
			if(!index.hint())
			{
				index.top(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, fl, fr);
				}

				index.left_levelzero(dim);
			}
		}

		index.get(dim, current_level, current_index);

		if (current_level > 0)
		{
			double fm = fml + fmr;

			double alpha_value = source[seq];

			double c = ((1/pow(2.0, static_cast<int>(current_level))) * static_cast<double>(current_index));

			// transposed operations:
			result[seq] = fm;

			fl = (fm/2.0) + (alpha_value*c) + fl;
			fr = (fm/2.0) - (alpha_value*c) + fr;
		}
		else
		{
			size_t seq_left;
			size_t seq_right;

			// left boundary
			seq_left = index.seq();

			// right boundary
			index.right_levelzero(dim);
			seq_right = index.seq();

			// up
			//////////////////////////////////////
			result[seq_left] = fl;
			result[seq_right] = fr;

			result[seq_left] -= 1.0/3.0*source[seq_right];

			index.left_levelzero(dim);
		}
	}
};

} // namespace detail

} // namespace sg

#endif /* SQXDPHIDPHIUPLINEARBOUNDARYUSCALED_HPP */
