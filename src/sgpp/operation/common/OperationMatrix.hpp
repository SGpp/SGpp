/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMATRIX_HPP
#define OPERATIONMATRIX_HPP

#include "data/DataVector.hpp"

namespace sg
{

/**
 * Multiplication with a Laplace matrix on sparse grids
 */
class OperationMatrix
{
public:
	/**
	 * Constructor
	 */
	OperationMatrix() {}

	/**
	 * Destructor
	 */
	virtual ~OperationMatrix() {}

	/**
	 * starts the Multiplication with the Laplace matrix
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 * @param result DataVector into which the result of the Laplace operation is stored
	 */
	virtual void mult(DataVector& alpha, DataVector& result) = 0;
};

}

#endif /* OPERATIONMATRIX_HPP */
