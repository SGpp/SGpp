/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
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
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef LAPLACEUPGRADIENTMODLINEAR_HPP
#define LAPLACEUPGRADIENTMODLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

namespace sg
{

namespace detail
{

class LaplaceUpGradientModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;

public:
	LaplaceUpGradientModLinear(GridStorage* storage) : storage(storage)
	{
	}

	~LaplaceUpGradientModLinear()
	{
	}

	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		double f = 0.0;
		rec(source, result, index, dim, f);
	}

protected:
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& f)
	{
		size_t seq = index.seq();

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double alpha_value = source[seq];
		double ht = pow(2.0, l);

		if(l == 1)
		{
			f = 0.0;
			if(!index.hint(dim))
			{
				index.left_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}

				f = 0.0;
				index.step_right(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}

			result[seq] = 0.0;
		}
		// left boundary
		else if(i == 1)
		{
			f = 0.0;
			if(!index.hint(dim))
			{
				index.left_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}

			result[seq] = ht * f;

			f += 2.0 * alpha_value;
		}
		// right boundary
		else if(i == (1 << l)-1)
		{
			f = 0.0;
			if(!index.hint(dim))
			{
				index.right_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}

			result[seq] = ht * f;

			f += 2.0 * alpha_value;
		}
	}

};


}

}

#endif /* LAPLACEUPGRADIENTMODLINEAR_HPP */
