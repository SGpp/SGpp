/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTEST_HPP
#define OPERATIONTEST_HPP

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

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
	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes) = 0;

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
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers) = 0;
};

}

#endif /* OPERATIONTEST_HPP */
