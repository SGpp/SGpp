/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIXSPAVXIDENTITY_HPP
#define DMSYSTEMMATRIXSPAVXIDENTITY_HPP

#include "data/DataVectorSP.hpp"
#include "grid/Grid.hpp"
#include "operation/datadriven/OperationBVectorizedSP.hpp"
#include "operation/common/OperationMatrixSP.hpp"

namespace sg
{

/**
 * Class that implements the virtual class OperationMatrixSP for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For Operation B an instance of OperationBVectorizedSPAVX is used!
 */
class DMSystemMatrixSPAVXIdentity : public OperationMatrixSP
{
private:
	/// the lambda, the regularisation parameter
	float lamb;
	/// OperationB for calculating the data matrix
	OperationBVectorizedSP* B;
	/// Pointer to the data vector
	DataMatrixSP* data;
	/// Number of orignal training instances
	size_t numTrainingInstances;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to DataMatrix that contains the training data
	 * @param C the regression functional
	 * @param lambda the lambda, the regression parameter
	 */
	DMSystemMatrixSPAVXIdentity(Grid& SparseGrid, DataMatrixSP& trainData, float lambda);

	/**
	 * Std-Destructor
	 */
	virtual ~DMSystemMatrixSPAVXIdentity();

	virtual void mult(DataVectorSP& alpha, DataVectorSP& result);

	/**
	 * Generates the right hand side of the classification equation
	 *
	 * @param classes the class information of the training data
	 * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
	 */
	void generateb(DataVectorSP& classes, DataVectorSP& b);

	/**
	 * rebuilds the DataMatrix for Level and Index
	 */
	void rebuildLevelAndIndex();
};

}

#endif /* DMSYSTEMMATRIXSPAVXIDENTITY_HPP */
