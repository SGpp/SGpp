/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
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

#ifndef LAPLACEUPMODLINEAR_HPP
#define LAPLACEUPMODLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

namespace sg
{

namespace detail
{

/**
 * Implements the up Method needed for the Laplace operator on mod linear grids
 */
class LaplaceUpModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	/// Pointer to the GridStorage object
	GridStorage* storage;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	LaplaceUpModLinear(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~LaplaceUpModLinear()
	{
	}

	/**
	 * This operations performs the calculation of up in the direction of dimension <i>dim</i>
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the up operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
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
	 * @param fl function value on the left boundary
	 * @param fr function value on the right boundary
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr)
	{
		size_t seq = index.seq();

		double alpha_value = source[seq];

		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double h = 1/pow(2.0, l);

		double fml = 0.0;
		double fmr = 0.0;

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

		double fm = fml + fmr;

		// level 1, constant function
		if(l == 1)
		{
			result[seq] = fl + fm + fr;

			fl += fm/2.0 + alpha_value;
			fr += fm/2.0 + alpha_value;
		}
		// left boundary
		else if(i == 1)
		{
			result[seq] = 2.0 * fl + fm;

			fl += fm/2.0 + 4.0/3.0*h*alpha_value;
			fr += fm/2.0 + 2.0/3.0*h*alpha_value;
		}
		// right boundary
		else if(static_cast<int>(i) == static_cast<int>((1 << l)-1))
		{
			result[seq] = 2.0 * fr + fm;

			fl += fm/2.0 + 2.0/3.0*h*alpha_value;
			fr += fm/2.0 + 4.0/3.0*h*alpha_value;
		}
		// inner functions
		else
		{
			result[seq] = fm;

			fl += fm/2.0 + h/2.0*alpha_value;
			fr += fm/2.0 + h/2.0*alpha_value;
		}
	}

};


}

}

#endif /* LAPLACEUPMODLINEAR_HPP */
