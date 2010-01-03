/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007-20009 Dirk Pflueger (dirk.pflueger@in.tum.de)          */
/* Copyright (C) 2007 Jörg Blank (blankj@in.tum.de)                          */
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

#ifndef LAPLACEUPLINEAR_HPP
#define LAPLACEUPLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * up-operation in dimension dim. for use with sweep, linear grids without boundaries
 */
class LaplaceUpLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// pointer to the grid's grid storage object
	GridStorage* storage;

public:
	/**
	 * Constructor
	 *
	 * @param storage pointer to the grid's grid storage object
	 */
	LaplaceUpLinear(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~LaplaceUpLinear()
	{
	}

	/**
	 * operator called by sweep during steping through the grid, start the calculation of Up
	 *
	 * @param source DataVector that contains the coefficients of the ansatzfunction
	 * @param result DataVector in which the result of the operation is stored
	 * @param index reference to a griditerator object that is used navigate through the grid
	 * @param dim the dimension in which the operation is executed
	 */
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		// provide memory for references
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

		{
			GridStorage::index_type::level_type l;
			GridStorage::index_type::index_type i;

			index.get(dim, l, i);

			double fm = fml + fmr;

			double alpha_value = source[seq];
			double h = 1/pow(2.0, static_cast<int>(l));

			// transposed operations:
			result[seq] = fm;

			fl = fm/2.0 + alpha_value*h/2.0 + fl;
			fr = fm/2.0 + alpha_value*h/2.0 + fr;
		}
	}
};

} // namespace detail

} // namespace sg

#endif /* LAPLACEUPLINEAR_HPP */
