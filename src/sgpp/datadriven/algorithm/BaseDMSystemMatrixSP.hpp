/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BASEDMSYSTEMMATRIXSP_HPP
#define BASEDMSYSTEMMATRIXSP_HPP

#include "base/datatypes/DataVectorSP.hpp"
#include "base/datatypes/DataMatrixSP.hpp"
#include "base/operation/OperationMatrixSP.hpp"

namespace sg
{
namespace datadriven
{

/**
 * Abstract class that defines the virtual class sg::base::OperationMatrix for
 * classification and regression problems (single precision version)
 */
class BaseDMSystemMatrixSP : public sg::base::OperationMatrixSP
{
private:
	/// the dataset
	sg::base::DataMatrixSP* data_;
	/// the lambda, the regularisation parameter
	double lambda_;

public:
	/**
	 * Std-Constructor
	 *
	 * @param lambda the lambda, the regression parameter
	 */
	BaseDMSystemMatrixSP(sg::base::DataMatrixSP& trainData, double lambda);

	/**
	 * Std-Destructor
	 */
	virtual ~BaseDMSystemMatrixSP();

	virtual void mult(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result) = 0;

	/**
	 * Generates the right hand side of the classification equation
	 *
	 * @param classes the class information of the training data
	 * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
	 */
	virtual void generateb(sg::base::DataVectorSP& classes, sg::base::DataVectorSP& b) = 0;
};

}
}

#endif /* BASEDMSYSTEMMATRIXSP_HPP */
