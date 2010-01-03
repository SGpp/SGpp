/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#ifndef OPERATIONTEST_HPP
#define OPERATIONTEST_HPP

#include "data/DataVector.hpp"

#ifdef WINDOWS
#pragma warning(disable: 4267)
#endif

namespace sg
{

/**
 * Operation the tests the function that is applied the current Sparse Grid at a given point
 *
 * @version $HEAD$
 */
class OperationTest
{
public:
	/**
	 * Constructor
	 */
	OperationTest() {}

	/**
	 * Destructor
	 */
	virtual ~OperationTest() {}

	/**
	 * Computes the classification accuracy on some test data.
	 *
	 * The function is evaluated at the given points. Tests on the classes {+1, -1}, cut-off at 0.
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 * @param data the coordinates of the evaluation points
	 * @param classes DataVector the holds the class information
	 */
	virtual double test(DataVector& alpha, DataVector& data, DataVector& classes) = 0;

	/**
	 * Computes the classification accuracy on some test data.
	 *
	 * The function is evaluated at the given points. Tests on the classes {+1, -1}, cut-off at 0.
	 *
	 * Also the number of the TP TN FP FN are determined
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 * @param data the coordinates of the evaluation points
	 * @param classes DataVector the holds the class information
	 * @param charaNumbers the number of true positives, true negatives, false positives, false negatives (Vector of length 4)
	 */
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataVector& data, DataVector& classes, DataVector& charaNumbers) = 0;
};

}

#endif /* OPERATIONTEST_HPP */
