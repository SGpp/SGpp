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
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef CONVERTPREWAVELETTOLINEAR_HPP
#define CONVERTPREWAVELETTOLINEAR_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"
#include <iostream>

namespace sg
{
namespace base
{



/**
 * Class that implements the transformation of a hierarchical prewavelet sparse grid to a
 * hierarchical linear sparse grid. Therefore the ()operator is implemented in order to use
 * the sweep algorithm for the grid traversal. Let the coefficients from the hat basis be
 * \f$ h_{l,i}\f$ and from the prewavelet basis \f$ u_{l,i} \f$. To calculate the surplusses,
 * temp values are needed:
 * \f[
 * (l,i)\neq G_{n}^{1}:t_{l,i}=-\frac{6}{10}u_{l,i\pm1}+t_{l+1,2i}
 * \f]
 * All temp values for levels greater than the maximal level of the grid are set to 0. The actual
 * transformation is calculated as follows:
 * \f{eqnarray*}{
 * h_{l,i}&=&u_{l,i}+\frac{1}{10}u_{l,i\pm2}+t_{l+1,2i}-\frac{1}{2}t_{l,i\pm1}\qquad\mbox{if \ensuremath{(l,i)} is an inner point}\\h_{l,i}&=&\frac{9}{10}u_{l,i}+\frac{1}{10}u_{l,i\pm2}+t_{l+1,2i}-\frac{1}{2}t_{l,i\pm1}\qquad\mbox{if \ensuremath{(l,i)} is at border}
 * \f}
 *
 *The picture depicts all needed variables in oder to perform the transformation:
 * \image html prewavelets_dehierarch.png "This picture shows all involved gridpoints (red crosses) and temp values (green circles) to calculate the new hierarchical coefficients (red arrows) and new temp values (green arrows)."
 */
class ConvertPrewaveletToLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	typedef GridStorage::index_type::level_type level_type;
	typedef GridStorage::index_type::index_type index_type;

	/// the grid object
	GridStorage* storage;



public:
	/**
	 * Constructor, must be bind to a grid
	 *
	 * @param storage the grid storage object of the the grid, on which the hierarchisation should be executed
	 */
	ConvertPrewaveletToLinear(GridStorage* storage) :
		storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~ConvertPrewaveletToLinear()
	{
	}

	/**
	 * Converts a given prewavelet base to a normal linear base.
	 */
	void operator()(DataVector& source, DataVector& result,
			grid_iterator& index, size_t dim)
	{
		level_type max_level = index.getGridDepth(dim);

		if (max_level == 1)
		{
			return;
		}

		level_type level = max_level;

		level_type init_level;
		index_type init_index;
		size_t _seq;
		size_t _seq_temp;
		double _val = 0.0;
		double* temp_current = 0;
		double* temp_old = 0;


		index.get(dim, init_level, init_index);

		for (; level > 1; --level)
		{
			if (level == max_level)
			{
				temp_current = new double[1 << level];
				temp_old = new double[1 << (level + 1)];
				for (int t = 0; t < (1 << (level + 1)); t++)
				{
					temp_old[t] = 0;
				}

				for (int t = 2; t < (1 << level); t = t + 2)
				{
					index.set(dim, level, t - 1);
					_seq = index.seq();
					_val = storage->end(_seq) ? 0.0 : source[_seq];
					temp_current[t] = -0.6 * _val;

					index.set(dim, level, t + 1);
					_seq = index.seq();
					_val = storage->end(_seq) ? 0.0 : source[_seq];
					temp_current[t] = temp_current[t] - 0.6 * _val;
				}
			}
			else
			{
				delete[] temp_old;
				temp_old = temp_current;
				temp_current = new double[(1 << level)];

				for (int t = 2; t < (1 << level); t = t + 2)
				{
					index.set(dim, level, t - 1);
					_seq = index.seq();
					_val = storage->end(_seq) ? 0.0 : source[_seq];
					temp_current[t] = -0.6 * _val + temp_old[t * 2];

					index.set(dim, level, t + 1);
					_seq = index.seq();
					_val = storage->end(_seq) ? 0.0 : source[_seq];
					temp_current[t] = temp_current[t] - 0.6 * _val;
				}

			}

			//Special treatment for first index
			index.set(dim, level, 1);
			_seq = index.seq();
			_val = storage->end(_seq) ? 0.0 : source[_seq];
			double current_value = _val;
			double left_value = 0;

			if (!storage->end(_seq))
				result[_seq] = 0.9 * current_value;

			index.set(dim, level, 3);
			_seq_temp = index.seq();
			_val = storage->end(_seq_temp) ? 0.0 : source[_seq_temp];

			if (!storage->end(_seq))
			{
				result[_seq] += 0.1 * _val;
				result[_seq] += temp_old[2] - 0.5 * temp_current[2];
			}

			for (int i = 3; i < (1 << level) - 2; i = i + 2)
			{
				index.set(dim, level, i);
				_seq = index.seq();

				index.set(dim, level, i + 2);
				_seq_temp = index.seq();

				left_value = current_value;
				_val = storage->end(_seq) ? 0.0 : source[_seq];
				current_value = _val;

				_val = storage->end(_seq_temp) ? 0.0 : source[_seq_temp];
				if (!storage->end(_seq))
				{
					result[_seq] = current_value + 0.1 * left_value + 0.1
							* _val;

					result[_seq] += temp_old[i * 2] - 0.5 * temp_current[i - 1]
							- 0.5 * temp_current[i + 1];
				}
			}

			//Special treatment for last index
			size_t last = (1 << level) - 1;
			index.set(dim, level, last);
			_seq = index.seq();
			_val = storage->end(_seq) ? 0.0 : source[_seq];
			if (!storage->end(_seq))
			{
				result[_seq] = 0.9 * _val + 0.1 * current_value;

				result[_seq] += temp_old[last * 2] - 0.5 * temp_current[last
						- 1];
			}

		}

		index.set(dim, init_level, init_index);
		_seq = index.seq();

		result[_seq] = source[_seq] + temp_current[2];

		delete[] temp_old;
		temp_old = 0;
		delete[] temp_current;
		temp_current = 0;

	}

};

 // namespace detail

} // namespace sg
}

#endif /* CONVERTPREWAVELETTOLINEAR_HPP */
