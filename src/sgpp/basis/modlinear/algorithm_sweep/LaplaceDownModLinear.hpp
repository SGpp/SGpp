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

#ifndef LAPLACEDOWNMODLINEAR_HPP
#define LAPLACEDOWNMODLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

namespace sg
{

namespace detail
{

/**
 * Implements the down Method needed for the Laplace operator on mod linear grids
 */
class LaplaceDownModLinear
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
	LaplaceDownModLinear(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~LaplaceDownModLinear()
	{
	}

	/**
	 * This operations performs the calculation of down in the direction of dimension <i>dim</i>
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the up operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		rec(source, result, index, dim, 0.0, 0.0);
	}

protected:

	/**
	 * recursive function for the calculation of Down
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

		double h = 1/pow(2.0, l);
		double fm;

		// level 1, constant function
		if(l == 1)
		{
			//integration
			result[seq] = 0.0 + alpha_value;

			//dehierarchisation
			fm = (fl + fr) / 2.0 + alpha_value;

			//boundary value
			fl += alpha_value;
			fr += alpha_value;
		}
		// left boundary
		else if(i == 1)
		{
			//integration
			result[seq] = 2.0/3.0 * h * (2.0*fl + fr)
                        + 8.0/3.0 * h * alpha_value;

            //dehierarchisation
            fm = (fl + fr) / 2.0 + alpha_value;

            //boundary value
            fl += 2.0 * alpha_value;
		}
		// right boundary
		else if(i == (1 << l)-1)
		{
			//integration
			result[seq] = 2.0/3.0 * h * (fl + 2.0*fr)
                        + 8.0/3.0 * h * alpha_value;

            //dehierarchisation
            fm = (fl + fr) / 2.0 + alpha_value;

            //boundary value
            fr += 2.0 * alpha_value;
		}
		// inner functions
		else
		{
			//integration
			result[seq] = h * (fl + fr) / 2.0
                       + 2.0/3.0 * h * alpha_value;

            //dehierarchisation
            fm = (fl + fr) / 2.0 + alpha_value;

			//boundary value

		}

		if(!index.hint(dim))
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


};


}

}

#endif /* LAPLACEDOWNMODLINEAR_HPP */
