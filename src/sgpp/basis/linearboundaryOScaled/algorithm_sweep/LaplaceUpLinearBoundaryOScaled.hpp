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

#ifndef LAPLACEUPLINEARBOUNDARYOSCALED_HPP
#define LAPLACEUPLINEARBOUNDARYOSCALED_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

namespace sg
{

namespace detail
{

/**
 * up-operation in dimension dim. for use with sweep
 */
class LaplaceUpLinearBoundaryOScaled
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;

public:
	LaplaceUpLinearBoundaryOScaled(GridStorage* storage) : storage(storage)
	{
	}

	~LaplaceUpLinearBoundaryOScaled()
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

		{
			GridStorage::index_type::level_type l;
			GridStorage::index_type::index_type i;

			index.get(dim, l, i);

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

#endif /* LAPLACEUPLINEARBOUNDARYOSCALED_HPP */
