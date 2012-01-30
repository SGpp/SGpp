/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BASEDMSYSTEMMATRIX_HPP
#define BASEDMSYSTEMMATRIX_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/operation/OperationMatrix.hpp"

namespace sg
{
namespace datadriven
{

class BaseDMSystemMatrix : public sg::base::OperationMatrix
{
private:
	/// the dataset
	sg::base::DataMatrix* data_;
	/// the lambda, the regularisation parameter
	double lambda_;

public:
	/**
	 * Std-Constructor
	 *
	 * @param lambda the lambda, the regression parameter
	 */
	BaseDMSystemMatrix(sg::base::DataMatrix& trainData, double lambda);

	/**
	 * Std-Destructor
	 */
	virtual ~Base::DMSystemMatrix();

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

	/**
	 * Generates the right hand side of the classification equation
	 *
	 * @param classes the class information of the training data
	 * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
	 */
	virtual void generateb(sg::base::DataVector& classes, sg::base::DataVector& b) = 0;
};

}
}

#endif /* BASEDMSYSTEMMATRIX_HPP */
