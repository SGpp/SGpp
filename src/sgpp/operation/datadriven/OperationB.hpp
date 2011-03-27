/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONB_HPP
#define OPERATIONB_HPP

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

namespace sg
{

/**
 * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
 * 
 * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data points, then B is a @f$N\times m@f$ matrix, with
 * @f[ (B)_{i,j} = \varphi_i(x_j). @f]
 */
class OperationB
{
public:
	/**
	 * Constructor
	 */
	OperationB() {}

	/**
	 * Destructor
	 */
	virtual ~OperationB() {}

	/**
	 * Multiplication of @f$B@f$ with vector @f$\alpha@f$
	 *
	 * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
	 * @param data vector, providing the data points x row-wise
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void mult(DataVector& alpha, DataMatrix& data, DataVector& result) = 0;

	/**
	 * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
	 *
	 * @param alpha vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
	 * @param data vector, providing the data points x row-wise
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result) = 0;
};

}

#endif /* OPERATIONB_HPP */
