/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIXAVXIDENTITY_HPP
#define DMSYSTEMMATRIXAVXIDENTITY_HPP

#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "operation/datadriven/OperationBVectorized.hpp"
#include "operation/common/OperationMatrix.hpp"

namespace sg
{

/**
 * Class that implements the virtual class OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For Operation B an instance of OperationBVectorizedAVX is used!
 */
class DMSystemMatrixAVXIdentity : public OperationMatrix
{
private:
	/// the lambda, the regularisation parameter
	double lamb;
	/// OperationB for calculating the data matrix
	OperationBVectorized* B;
	/// Pointer to the data vector
	DataMatrix* data;
	/// Number of orignal training instances
	size_t numTrainingInstances;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to DataVector that contains the training data
	 * @param C the regression functional
	 * @param lambda the lambda, the regression parameter
	 */
	DMSystemMatrixAVXIdentity(Grid& SparseGrid, DataMatrix& trainData, double lambda);

	/**
	 * Std-Destructor
	 */
	virtual ~DMSystemMatrixAVXIdentity();

	virtual void mult(DataVector& alpha, DataVector& result);

	/**
	 * Generates the right hand side of the classification equation
	 *
	 * @param classes the class information of the training data
	 * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
	 */
	void generateb(DataVector& classes, DataVector& b);

	/**
	 * rebuilds the DataMatrix for Level and Index
	 */
	void rebuildLevelAndIndex();
};

}

#endif /* DMSYSTEMMATRIXAVXIDENTITY_HPP */
