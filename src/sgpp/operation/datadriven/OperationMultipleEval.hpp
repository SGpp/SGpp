/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVAL_HPP
#define OPERATIONMULTIPLEEVAL_HPP

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "operation/common/OperationMatrix.hpp"

namespace sg
{

/**
 * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
 *
 * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data points
 */
class OperationMutlipleEval
{
protected:
	/// Pointer to the DataSet that should be evaluated on the grid
	DataMatrix* dataset_;

public:
	/**
	 * Constructor
	 */
	OperationMutlipleEval(DataMatrix& DataSet) {}

	/**
	 * Destructor
	 */
	virtual ~OperationMutlipleEval() {}

	/**
	 * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
	 *
	 * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void mult(DataVector& alpha, DataVector& result) = 0;

	/**
	 * Multiplication of @f$B@f$ with vector @f$\alpha@f$
	 *
	 * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void multTranspose(DataVector& source, DataVector& result) = 0;
};

}

#endif /* OPERATIONMULTIPLEEVAL_HPP */
