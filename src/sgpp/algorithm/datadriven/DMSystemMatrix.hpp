/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIX_HPP
#define DMSYSTEMMATRIX_HPP

#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "operation/datadriven/OperationB.hpp"
#include "operation/common/OperationMatrix.hpp"

namespace sg
{

/**
 * Class that implements the virtual class OperationMatrix for the
 * application of classification for the Systemmatrix
 */
class DMSystemMatrix : public OperationMatrix
{
private:
	/// the lambda, the regularisation parameter
	double lamb;
	/// OperationMatrix, the regularisation mehtod
	OperationMatrix* C;
	/// OperationB for calculating the data matrix
	OperationB* B;
	/// Pointer to the data vector
	DataMatrix* data;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to DataVector that contains the training data
	 * @param C the regression functional
	 * @param lambda the lambda, the regression parameter
	 */
	DMSystemMatrix(Grid& SparseGrid, DataMatrix& trainData, OperationMatrix& C, double lambda);

	/**
	 * Std-Destructor
	 */
	virtual ~DMSystemMatrix();

	virtual void mult(DataVector& alpha, DataVector& result);

	/**
	 * Generates the right hand side of the classification equation
	 *
	 * @param classes the class information of the training data
	 * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
	 */
	void generateb(DataVector& classes, DataVector& b);
};

}

#endif /* DMSYSTEMMATRIX_HPP */
