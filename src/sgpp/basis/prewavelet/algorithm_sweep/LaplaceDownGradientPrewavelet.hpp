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

#ifndef LAPLACEDOWNGRADIENTPREWAVELET_HPP
#define LAPLACEDOWNGRADIENTPREWAVELET_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{
namespace pde
{



/**
 * Implements the downGradient Method needed for the Laplace operator on prewavelet grids.
 * The calculation is done iterative and utilizes the following temp variables:
 * \f{eqnarray*}{
 * t_{k,j}&=&\frac{1}{h_{k}}\left(\frac{32}{10}u_{k,j}+\frac{8}{10}u_{k,j\pm2}\right) \quad(k,j)\in G_{n}^{1}\\ t_{k,j}&=&\frac{1}{h_{k}}\left(\frac{24}{10}u_{k,j}+\frac{8}{10}u_{k,j\pm2}\right) \quad(k,j)\in G_{n}^{1}\mbox{ and next to border}\\ t_{k,j}&=&t_{k-1,\frac{j}{2}}-\frac{1}{h_{k}}\left(\frac{23}{10}u_{k,j\pm1}+\frac{1}{10}u_{k,j\pm3}\right) \quad(k,j)\notin G_{n}^{1}\\ t_{k,j}&=&t_{k-1,\frac{j}{2}}-\frac{1}{h_{k}}\left(\frac{22}{10}u_{k,j\pm1}+\frac{1}{10}u_{k,j\pm3}\right) \quad(k,j)\notin G_{n}^{1}\mbox{ and \ensuremath{(k,j\pm1)}next to border}
 * \f}
 *
 * The calculation of the results is as follows:
 * \f{eqnarray*}{
 r_{k,j}&=&-\frac{6}{10}t_{k-1,\frac{j\pm1}{2}}\\
 &&+\frac{1}{h_{k}}u_{k,j}\cdot\left[\frac{612}{100}\mbox{ if \ensuremath{(k,j)} inner point, or }
 \frac{356}{100}\mbox{ if \ensuremath{(k,j)} border point}\right]\\
 &&+\frac{1}{h_{k}}u_{k,j\pm2}\cdot\left[\frac{256}{100} \mbox{ if \ensuremath{(k,j)} and \ensuremath{(k,j\pm2)} inner point, or }
 \frac{242}{100} \mbox{ if \ensuremath{(k,j)} or \ensuremath{(k,j\pm2)} border point, or }
 \frac{228}{100} \mbox{ if \ensuremath{(k,j)} border point}\right]\\
 &&+\frac{1}{h_{k}}u_{k,j\pm4}\cdot\frac{14}{100}\\r_{1,1}&=&4u_{1,1}
 * \f}
 * Please note, that all values of gridpoints outside of the sparse grid are treated as 0. The following
 * picture depicts all involved grid points and temp values in order to calculate a specific point:
 * \image html prewavelets_down.png "All involved gridpoint for the up algorithm (red) and temp points between grid points (green). The gray line indicates the support of the prewavelet."
 */

class LaplaceDownGradientPrewavelet
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
	LaplaceDownGradientPrewavelet(sg::base::GridStorage* storage) :
		storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~LaplaceDownGradientPrewavelet()
	{
	}

	/**
	 * This operations performs the calculation of downGradient in the direction of dimension <i>dim</i>
	 *
	 * @param source sg::base::DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result sg::base::DataVector that contains the result of the down operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	void operator()(sg::base::DataVector& source, sg::base::DataVector& result,
			grid_iterator& index, size_t dim)

	{

		size_t _seql2, _seql1, _seq, _seqr1, _seqr2 = 0;
		double _vall2, _vall1, _val, _valr1, _valr2 = 0;

		sg::base::GridStorage::index_type::level_type l, l_old, max_level;
		sg::base::GridStorage::index_type::index_type i, i_old;

		double h;

		index.get(dim, l_old, i_old); //save iterator and restore afterwards
		max_level = index.getGridDepth(dim); //How many levels to decent?

		//We need the temp_values of the level above
		double* temp_current = new double[1];
		double* temp_old;

		//Level 1-----------------------------------------------------------
		_seq = index.seq();
		temp_current[0] = 4 * source[_seq];
		result[_seq] = 4 * source[_seq];

		//Is there a second level?
		if (max_level == 1)
		{
			//No, delete temp and return
			delete[] temp_current;
			index.set(dim, l_old, i_old);
			return;
		}

		//Level 2-----------------------------------------------------------
		l = 2;
		h = 1 << l;
		temp_old = temp_current;
		temp_current = new double[3];
		index.left_child(dim);
		_seql1 = index.seq();
		index.step_right(dim);
		_seqr1 = index.seq();

		//Calculate new temp variables
		temp_current[0] = 2.4 * h * source[_seql1] + 0.8 * h * source[_seqr1];
		temp_current[1] = temp_old[0] - 2.2 * h * (source[_seql1]
				+ source[_seqr1]);
		temp_current[2] = 2.4 * h * source[_seqr1] + 0.8 * h * source[_seql1];

		//Calculate results
		_vall1 = source[_seql1];//We save the variables to enable work in sito.
		_valr1 = source[_seqr1];
		result[_seql1] = -0.6 * temp_old[0] + 3.56 * h * _vall1 + 2.28 * h
				* _valr1;
		result[_seqr1] = -0.6 * temp_old[0] + 3.56 * h * _valr1 + 2.28 * h
				* _vall1;

		//All other levels--------------------------------------------------
		for (l = 3; l <= max_level; l++)
		{
			h = 1 << l;
			index.set(dim, l, 1);

			//Calculate new temp-Variables ################################################
			delete[] temp_old;
			temp_old = temp_current;
			temp_current = new double[(1 << l) - 1];

			if (l != max_level) //If we are in the last iteration, we don't need any more temp-variables, thus skip the calculation
			{
				//Ramp up
				_seql2 = index.seq();
				_vall2 = storage->end(_seql2) ? 0.0 : source[_seql2];

				index.step_right(dim);
				_seql1 = index.seq();
				_vall1 = storage->end(_seql1) ? 0.0 : source[_seql1];

				index.step_right(dim);
				_seqr1 = index.seq();
				_valr1 = storage->end(_seqr1) ? 0.0 : source[_seqr1];

				index.step_right(dim);
				_seqr2 = index.seq();
				_valr2 = storage->end(_seqr2) ? 0.0 : source[_seqr2];

				temp_current[0] = 2.4 * h * _vall2// sg::base::Grid-point
						+ 0.8 * h * _vall1; //right neighbor

				temp_current[1] = temp_old[0]//
						- 2.2 * h * _vall2 // left neighbor
						- 2.3 * h * _vall1 // right neighbor
						- 0.1 * h * _valr1; //right-right neighbor

				for (i = 2; (int)i < (1 << l) - 3; i++)
				{
					if (i % 2 == 0) //On sg::base::Grid-Points
					{
						if ((int)i == (1 << l) - 4)
						{
							temp_current[i] = 3.2 * h * _valr1 //current grid point
									+ 0.8 * h * (_vall1 + _valr2); //neighbors left and right
						}
						else
						{
							temp_current[i] = 3.2 * h * _vall1 //current grid point
									+ 0.8 * h * (_vall2 + _valr1); //neighbors left and right
						}
					}
					else //between grid points
					{
						temp_current[i] = temp_old[(i - 1) / 2] //temp-value from above
								- 2.3 * h * (_vall1 + _valr1)//left and right neighbor
								- 0.1 * h * (_vall2 + _valr2);//left-left and right-right neighbor

						//forward indoces except in the last iteration
						if ((int)i != (1 << l) - 5)
						{
							_vall2 = _vall1;
							_vall1 = _valr1;
							_valr1 = _valr2;
							index.step_right(dim);
							_seq = index.seq();
							_valr2 = storage->end(_seq) ? 0.0 : source[_seq];
						}

					}
				}

				//Two temp-values left
				temp_current[(1 << l) - 3] = temp_old[((1 << l) - 4) / 2] //
						- 2.2 * h * _valr2 - 2.3 * h * _valr1 //left and right neighbor
						- 0.1 * h * _vall1; //left-left neighbor

				temp_current[(1 << l) - 2] = 2.4 * h * _valr2 + 0.8 * h
						* _valr1; //On sg::base::Grid-point


			}
			// End of temp-value calculation ################################################

			// Calculate result ################################################
			//Ramp-up
			index.set(dim, l, 1);
			_seql2 = index.seq();
			_vall2 = storage->end(_seql2) ? 0.0 : source[_seql2];

			index.step_right(dim);
			_seql1 = index.seq();
			_vall1 = storage->end(_seql1) ? 0.0 : source[_seql1];

			index.step_right(dim);
			_seq = index.seq();
			_val = storage->end(_seq) ? 0.0 : source[_seq];

			index.step_right(dim);
			_seqr1 = index.seq();
			_valr1 = storage->end(_seqr1) ? 0.0 : source[_seqr1];

			if (l != 3)
			{
				index.step_right(dim);
				_seqr2 = index.seq();
				_valr2 = storage->end(_seqr2) ? 0.0 : source[_seqr2];
			}

			if (!storage->end(_seql2))
				result[_seql2] = -0.6 * temp_old[0] + 3.56 * h * _vall2 //current point
						+ 2.42 * h * _vall1 //right neighbor
						+ 0.14 * h * _val; //right-right neighbor

			if (!storage->end(_seql1))
				result[_seql1] = -0.6 * (temp_old[0] + temp_old[1]) //
						+ 6.12 * h * _vall1 //current point
						+ 2.42 * h * _vall2 //left neighbor
						+ 2.56 * h * _val //right neighbor
						+ 0.14 * h * _valr1; //right-right neighbor

			if (l == 3)
			{
				/* level = 3, that means we have only 4 gridpoints, thus we have to
				 * shift the indices so that ramp-down work correctly. Otherwise
				 * proceed with the inner-points
				 */
				_valr2 = _valr1;
				_seqr2 = _seqr1;

				_valr1 = _val;
				_seqr1 = _seq;

				_val = _vall1;
				_seq = _seql1;

				_vall1 = _vall2;
				_seql1 = _seql2;
			}
			else
			{

				for (i = 5; (int)i < (1 << l) - 4; i += 2)
				{
					if (!storage->end(_seq))
						result[_seq] = -0.6 * (temp_old[(i - 3) / 2]
								+ temp_old[(i - 1) / 2]) //
								+ 6.12 * h * _val //current point
								+ 2.56 * h * (_vall1 + _valr1) //left and right neighbor
								+ 0.14 * h * (_vall2 + _valr2); //left-left and right-right neighbor

					if ((int)i != (1 << l) - 5)
					{
						//forward indices, except during the last iteration
						_vall2 = _vall1;
						_vall1 = _val;
						_val = _valr1;
						_seq = _seqr1;
						_valr1 = _valr2;
						_seqr1 = _seqr2;
						index.step_right(dim);
						_seqr2 = index.seq();
						_valr2 = storage->end(_seqr2) ? 0.0 : source[_seqr2];
					}
				}
			}

			//Again, two points left
			if (!storage->end(_seqr1))
				result[_seqr1] = -0.6 * (temp_old[(1 << (l - 1)) - 3]
						+ temp_old[(1 << (l - 1)) - 2]) //
						+ 6.12 * h * _valr1 //current point
						+ 2.42 * h * _valr2 //right neighbor
						+ 2.56 * h * _val //left neighbor
						+ 0.14 * h * _vall1; //left-left neighbor

			if (!storage->end(_seqr2))
				result[_seqr2] = -0.6 * temp_old[(1 << (l - 1)) - 2]//
						+ 3.56 * h * _valr2 //current point
						+ 2.42 * h * _valr1 //left neighbor
						+ 0.14 * h * _val; //left-left neighbor

		}

		index.set(dim, l_old, i_old);
		delete[] temp_current;
		temp_current = 0;
		delete[] temp_old;
		temp_old = 0;
	}

};



}
}

#endif /* LAPLACEDOWNGRADIENTPREWAVELET_HPP */
