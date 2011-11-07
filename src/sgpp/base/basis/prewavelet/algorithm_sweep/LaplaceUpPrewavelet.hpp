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

#ifndef LAPLACEUPPREWAVELET_HPP
#define LAPLACEUPPREWAVELET_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg
{
namespace pde
{



/**
 * Implements the down Method needed for the Laplace operator on prewavelet grids.
 * Here, the prewavlets have an advantage over the normal linear base due to the semi
 * othogonality. That means, for the calculation, only the gridpoint in the same level have
 * to be taken into account. The matrix entries are the following:
 * \f{eqnarray*}{
 * m_{(1,1)(1,1)}&=&\frac{1}{3}\\m_{(2,1)(2,3)}&=&\frac{1}{25}\\m_{(l,1)(l,1)}&=&m_{(l,2^{l}-1)(l,2^{l}-1)}=\frac{44}{75}h_{l}\\m_{(l,1)(l,3)}&=&m_{(l,2^{l}-1)(l,2^{l}-3)}=\frac{11}{75}h_{l}\\m_{(l,i)(l,i)}&=&\frac{18}{25}h_{l}\\m_{(l,i)(l,i\pm2)}&=&\frac{2}{15}h_{l}\\m_{(l,i)(l,i\pm4)}&=&-\frac{11}{75}h_{l}
 * \f}
 * With that, the calculation of the entire matrix is completed, thus no additional up-part is needed.
 */
class LaplaceUpPrewavelet
{
protected:
	typedef sg::base::GridStorage::grid_iterator grid_iterator;
	/// Pointer to sg::base::GridStorage object
	sg::base::GridStorage* storage;

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's sg::base::GridStorage object
	 */
	LaplaceUpPrewavelet(sg::base::GridStorage* storage) :
		storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~LaplaceUpPrewavelet()
	{
	}

	/**
	 * This operations performs the calculation of down in the direction of dimension <i>dim</i>
	 *
	 * @param source sg::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result sg::base::DataVector that contains the result of the up operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(sg::base::DataVector& source, sg::base::DataVector& result,
			grid_iterator& index, size_t dim)
	{

		size_t seq = index.seq();
		sg::base::GridStorage::index_type::level_type l;
		sg::base::GridStorage::index_type::index_type i;
		sg::base::GridStorage::index_type::level_type l_old;
		sg::base::GridStorage::index_type::index_type i_old;
		sg::base::GridStorage::index_type::index_type last_index;
		size_t _seq;
		size_t _seql1;
		size_t _seql2;
		size_t _seqr1;
		size_t _seqr2;
		double _val, _vall1, _vall2, _valr1, _valr2;
		double h;
		bool hasChilds = false;
		//GridStorage::index_type::level_type max_level = getGridDepth(index, dim);

		index.get(dim, l, i);
		index.get(dim, l_old, i_old);

		//Level 1

		result[seq] = 1.0 / 3.0 * source[seq];

		if (!index.hint_left(dim) && !index.hint_right(dim))
		{
			return;
		}

		//Level 2
		l = 2;
		h = 1 << l;
		index.set(dim, 2, 1);
		if (!hasChilds && (index.hint_left(dim) || index.hint_right(dim)))
			hasChilds = true;
		_seql1 = index.seq();
		_vall1 = storage->end(_seql1) ? 0.0 : source[_seql1];

		index.set(dim, 2, 3);
		if (!hasChilds && (index.hint_left(dim) || index.hint_right(dim)))
			hasChilds = true;
		_seqr1 = index.seq();
		_valr1 = storage->end(_seqr1) ? 0.0 : source[_seqr1];

		if (!storage->end(_seql1))
			result[_seql1] = 11.0 / 75.0 * _vall1 + 1.0 / 25.0 * _valr1;

		if (!storage->end(_seqr1))
			result[_seqr1] = 11.0 / 75.0 * _valr1 + 1.0 / 25.0 * _vall1;
		//
		//		result[_seql1] = 11.0 / 75.0 * source[_seql1];
		//		result[_seqr1] = 11.0 / 75.0 * source[_seqr1];


		if (!hasChilds)
		{
			index.set(dim, l_old, i_old);
			return;
		}

		while (true)
		{
			l++;
			hasChilds = false;
			h = 1.0 / (1 << l);

			last_index = (1 << (l - 1)) - 1; //Number of Points in this level

			index.set(dim, l, 1);
			if (!hasChilds && (index.hint_left(dim) || index.hint_right(dim)))
				hasChilds = true;
			_seq = index.seq();
			_val = storage->end(_seq) ? 0.0 : source[_seq];

			index.set(dim, l, 3);
			if (!hasChilds && (index.hint_left(dim) || index.hint_right(dim)))
				hasChilds = true;
			_seqr1 = index.seq();
			_valr1 = storage->end(_seqr1) ? 0.0 : source[_seqr1];

			index.set(dim, l, 5);
			if (!hasChilds && (index.hint_left(dim) || index.hint_right(dim)))
				hasChilds = true;
			_seqr2 = index.seq();
			_valr2 = storage->end(_seqr2) ? 0.0 : source[_seqr2];

			if (!storage->end(_seq))
				result[_seq] = 44.0 / 75.0 * h * _val //
						+ 11.0 / 75.0 * h * _valr1//
						- 1.0 / 75.0 * h * _valr2;//

			_seql1 = _seq;
			_seq = _seqr1;
			_seqr1 = _seqr2;
			_vall1 = _val;
			_val = _valr1;
			_valr1 = _valr2;
			index.set(dim, l, 7);
			if (!hasChilds && (index.hint_left(dim) || index.hint_right(dim)))
				hasChilds = true;
			_seqr2 = index.seq();
			_valr2 = storage->end(_seqr2) ? 0.0 : source[_seqr2];

			if (!storage->end(_seq))
				result[_seq] = 11.0 / 75.0 * h * _vall1 //
						+ 18.0 / 25.0 * h * _val //
						+ 2.0 / 15.0 * h * _valr1//
						- 1.0 / 75.0 * h * _valr2;//

			//Main loop--------------------------------------

			for (i = 2; i < last_index - 1; i++)
			{

				_seql2 = _seql1;
				_seql1 = _seq;
				_seq = _seqr1;
				_seqr1 = _seqr2;
				_vall2 = _vall1;
				_vall1 = _val;
				_val = _valr1;
				_valr1 = _valr2;
				index.set(dim, l, i * 2 + 5);
				if (!hasChilds && (index.hint_left(dim)
						|| index.hint_right(dim)))
					hasChilds = true;
				_seqr2 = index.seq();
				_valr2 = storage->end(_seqr2) ? 0.0 : source[_seqr2];

				if (!storage->end(_seq))
					result[_seq] = -1.0 / 75.0 * h * _vall2 //
							+ 2.0 / 15.0 * h * _vall1 //
							+ 18.0 / 25.0 * h * _val //
							+ 2.0 / 15.0 * h * _valr1//
							- 1.0 / 75.0 * h * _valr2;//

			}

			//Main loop--------------------------------------

			_seql2 = _seql1;
			_seql1 = _seq;
			_seq = _seqr1;
			_seqr1 = _seqr2;
			_vall2 = _vall1;
			_vall1 = _val;
			_val = _valr1;
			_valr1 = _valr2;

			if (!storage->end(_seq))
			result[_seq] = 11.0 / 75.0 * h * _valr1 //
					+ 18.0 / 25.0 * h * _val //
					+ 2.0 / 15.0 * h * _vall1//
					- 1.0 / 75.0 * h * _vall2;//

			_seql2 = _seql1;
			_seql1 = _seq;
			_seq = _seqr1;
			_vall2 = _vall1;
			_vall1 = _val;
			_val = _valr1;

			if (!storage->end(_seq))
			result[_seq] = 44.0 / 75.0 * h * _val //
					+ 11.0 / 75.0 * h * _vall1//
					- 1.0 / 75.0 * h * _vall2;//

			if (!hasChilds)
			{
				index.set(dim, l_old, i_old);
				return;
			}

		}

	}

protected:

};



}
}

#endif /* LAPLACEUPPREWAVELET_HPP */
