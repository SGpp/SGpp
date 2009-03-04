/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef LAPLACEUPLINEARBOUNDARY_HPP
#define LAPLACEUPLINEARBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

namespace sg
{

namespace detail
{

/**
 * up-operation in dimension dim. for use with sweep
 */
class LaplaceUpLinearBoundary
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;

public:
	LaplaceUpLinearBoundary(GridStorage* storage) : storage(storage)
	{
	}

	~LaplaceUpLinearBoundary()
	{
	}

	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		// get boundary values
		double fl = 0.0;
		double fr = 0.0;

		rec(source, result, index, dim, fl, fr);
	}

protected:

	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr)
	{
		size_t seq = index.seq();

		fl = fr = 0.0;
		double fml = 0.0;
		double fmr = 0.0;

		GridStorage::index_type::level_type current_level;
		GridStorage::index_type::index_type current_index;

		index.get(dim, current_level, current_index);

		if(current_level == 0)
		{
			if(!index.hint(dim))
			{
				index.top(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, fl, fr);
				}
				index.left_levelzero(dim);
			}
		}
		else
		{
			if(!index.hint(dim))
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

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		if (l == 0)
		{
			double alpha_value = 0.0;
			double fm = fml + fmr;

			// handle left boundary
			alpha_value = source[seq];

			// transposed operations:
			result[seq] = fm/2.0 + alpha_value/2.0 + fl;

			// handle right boundary
			index.right_levelzero(dim);
			seq = index.seq();

			alpha_value = source[seq];

			// transposed operations:
			result[seq] = fm/2.0 + alpha_value/2.0 + fr;

			index.left_levelzero(dim);
		}
		else
		{
			double fm = fml + fmr;

			double alpha_value = source[seq];
			double h = 1/pow(2.0,l);

			// transposed operations:
			result[seq] = fm;

			fl = fm/2.0 + alpha_value*h/2.0 + fl;
			fr = fm/2.0 + alpha_value*h/2.0 + fr;
		}
	}

};

} // namespace detail

} // namespace sg

#endif /* LAPLACEUPLINEARBOUNDARY_HPP */
