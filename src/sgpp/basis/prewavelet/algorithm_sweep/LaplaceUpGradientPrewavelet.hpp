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

#ifndef LAPLACEUPGRADIENTPREWAVELET_HPP
#define LAPLACEUPGRADIENTPREWAVELET_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"
using namespace sg::base;

namespace sg
{
namespace pde
{



/**
 * Implements the upGradient Method needed for the Laplace operator on prewavelet grids. The calculation
 * is done iterative and utilizes the following temp variables:
 * \f[
 * t_{k,j}=-\frac{6}{10}u_{k,j\pm1}+t_{k+1,2j}\qquad(k,j)\notin G_{n}^{1}
 * \f]
 * The correct values are then calculated as follows:
 * \f{eqnarray*}{
 * r_{k,j}&=&\frac{1}{h_{k}}\left(2t_{k+1,2j}-t_{k+1,2(j\pm1)}\right)\\&&-\frac{6}{10}\frac{1}{h_{k}}\left(-t_{k+1,2(j\pm2)}+2t_{k+1,2(j\pm1)}-2t_{k+1,2j}\right)\\&&+\frac{1}{10}\frac{1}{h_{k}}\left(-t_{k+1,2(j\pm1)}+2t_{k+1,2(j\pm2)}-t_{k+1,2(j\pm2)}\right)
 * \f}
 * In case of borders:
 * \f{eqnarray*}{
 * r_{k,j}&=&\frac{9}{10}\frac{1}{h_{k}}\left(2t_{k+1,2j}-t_{k+1,2(j\pm1)}\right)\\&&-\frac{6}{10}\frac{1}{h_{k}}\left(-t_{k+1,2(j\pm2)}+2t_{k+1,2(j\pm1)}-t_{k+1,2j}\right)\\&&+\frac{1}{10}\frac{1}{h_{k}}\left(-t_{k+1,2(j\pm1)}+2t_{k+1,2(j\pm2)}-t_{k+1,2(j\pm2)}\right)
 * \f}
 * Please note, that all values of gridpoints outside of the sparse grid are treated as 0. The following
 * picture depicts all involved grid points and temp values in order to calculate a specific point:
 * \image html prewavelets_up.png "All involved gridpoint for the up algorithm (red) and temp points between grid points (green). The gray line indicates the support of the prewavelet."
 */
class LaplaceUpGradientPrewavelet
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
	LaplaceUpGradientPrewavelet(GridStorage* storage) :
		storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~LaplaceUpGradientPrewavelet()
	{
	}

	/**
	 * This operations performs the calculation of upGradient in the direction of dimension <i>dim</i>
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the down operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(DataVector& source, DataVector& result,
			grid_iterator& index, size_t dim)
	{
		GridStorage::index_type::level_type l = index.getGridDepth(dim);
		GridStorage::index_type::index_type i;

		GridStorage::index_type::level_type l_old;
		GridStorage::index_type::index_type i_old;

		index.get(dim, l_old, i_old);

		size_t _seq;
		size_t _seql;
		size_t _seqr;

		double _vall, _valr;

		double h;

		if (l == 1)
			return;

		index.set(dim, l, 1);
		_seqr = index.seq();
		_valr = storage->end(_seqr) ? 0.0 : source[_seqr];

		double* temp_current = new double[(1 << (l - 1)) - 1];
		for (i = 0; i < (unsigned int) (1 << (l - 1)) - 1; i++)
		{
			_seql = _seqr;
			_vall = _valr;
			index.set(dim, l, 2 * i + 3);
			_seqr = index.seq();
			_valr = storage->end(_seqr) ? 0.0 : source[_seqr];
			temp_current[i] = -0.6 * (_vall + _valr);
		}
		l--;

		for (; l > 2; l--)
		{
			i = 0;
			h = 1 << l;
			index.set(dim, l, i + 1);
			_seq = index.seq();
			if (!storage->end(_seq))
				result[_seq] = 0.9 * h * (2 * temp_current[i] - temp_current[i
						+ 1]) - 0.6 * h * (-temp_current[i + 2] + 2
						* temp_current[i + 1] - temp_current[i]) + 0.1 * h
						* (-temp_current[i + 1] + 2 * temp_current[i + 2]
								- temp_current[i + 3]);

			i = 2;
			index.set(dim, l, i + 1);
			_seq = index.seq();
			if (!storage->end(_seq))
				result[_seq] = h * (2 * temp_current[i] - temp_current[i + 1]
						- temp_current[i - 1]) + -0.6 * h * (-temp_current[i
						+ 2] + 2 * temp_current[i + 1] - 2 * temp_current[i])
						- 0.6 * h * (-temp_current[i - 2] + 2 * temp_current[i
								- 1]) + 0.1 * h * (-temp_current[i + 1] + 2
						* temp_current[i + 2] - temp_current[i + 3]) + 0.1 * h
						* (-temp_current[i - 1] + 2 * temp_current[i - 2]);

			for (i = 4; i < (unsigned int) (1 << l) - 4; i = i + 2)
			{
				index.set(dim, l, i + 1);
				_seq = index.seq();
				if (!storage->end(_seq))
					result[_seq] = h * (2 * temp_current[i] - temp_current[i
							+ 1] - temp_current[i - 1]) - 0.6 * h
							* (-temp_current[i + 2] + 2 * temp_current[i + 1]
									- 2 * temp_current[i]) - 0.6 * h
							* (-temp_current[i - 2] + 2 * temp_current[i - 1])
							+ 0.1 * h
									* (-temp_current[i + 1] + 2
											* temp_current[i + 2]
											- temp_current[i + 3]) + 0.1 * h
							* (-temp_current[i - 1] + 2 * temp_current[i - 2]
									- temp_current[i - 3]);
			}

			i = (1 << l) - 4;
			index.set(dim, l, i + 1);
			_seq = index.seq();
			if (!storage->end(_seq))
				result[_seq] = h * (2 * temp_current[i] - temp_current[i + 1]
						- temp_current[i - 1]) + -0.6 * h * (-temp_current[i
						+ 2] + 2 * temp_current[i + 1] - 2 * temp_current[i])
						- 0.6 * h * (-temp_current[i - 2] + 2 * temp_current[i
								- 1]) + 0.1 * h * (-temp_current[i + 1] + 2
						* temp_current[i + 2]) + 0.1 * h
						* (-temp_current[i - 1] + 2 * temp_current[i - 2]
								- temp_current[i - 3]);

			i = (1 << l) - 2;
			index.set(dim, l, i + 1);
			_seq = index.seq();
			if (!storage->end(_seq))
				result[_seq] = 0.9 * h * (2 * temp_current[i] - temp_current[i
						- 1]) - 0.6 * h * (-temp_current[i - 2] + 2
						* temp_current[i - 1] - temp_current[i]) + 0.1 * h
						* (-temp_current[i - 1] + 2 * temp_current[i - 2]
								- temp_current[i - 3]);

			index.set(dim, l, 1);
			_seqr = index.seq();
			_valr = storage->end(_seqr) ? 0.0 : source[_seqr];
			for (i = 0; i < (unsigned int) (1 << (l - 1)) - 1; i++)
			{

				_seql = _seqr;
				_vall = _valr;
				index.set(dim, l, 2 * i + 3);
				_seqr = index.seq();
				_valr = storage->end(_seqr) ? 0.0 : source[_seqr];
				temp_current[i] = -0.6 * (_vall + _valr) + temp_current[2 * i
						+ 1];
			}
		}

		if (l == 2)
		{
			h = 1 << l;

			index.set(dim, 2, 1);
			_seql = index.seq();
			if (!storage->end(_seql)){
				result[_seql] = 0.9 * h * (2 * temp_current[0]
						- temp_current[1]) - 0.6 * h * (-temp_current[2] + 2
						* temp_current[1] - temp_current[0]) + 0.1 * h
						* (-temp_current[1] + 2 * temp_current[2]);
				_vall = source[_seql];
			}else{
				_vall = 0.0;
			}
			index.set(dim, 2, 3);
			_seqr = index.seq();
			if (!storage->end(_seqr)){
				result[_seqr] = 0.9 * h * (2 * temp_current[2]
						- temp_current[1]) - 0.6 * h * (-temp_current[0] + 2
						* temp_current[1] - temp_current[2]) + 0.1 * h
						* (-temp_current[1] + 2 * temp_current[0]);
				_valr = source[_seqr];
			}else
				_valr = 0.0;

			temp_current[0] = -0.6 * (_vall + _valr)
					+ temp_current[1];

			l = 1;
		}

		if (l == 1)
		{
			h = 1 << l;
			index.set(dim, 1, 1);
			_seq = index.seq();
			result[_seq] = 2 * h * temp_current[0];
		}

		//I dont think thats necessary, but just to be sure!
		index.set(dim, l_old, i_old);
		delete[] temp_current;

	}

};



}
}

#endif /* LAPLACEUPGRADIENTMODLINEAR_HPP */
