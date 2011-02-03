/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALVECTORIZEDSP_HPP
#define OPERATIONMULTIPLEEVALVECTORIZEDSP_HPP

#include "data/DataVectorSP.hpp"
#include "data/DataMatrixSP.hpp"
#include "operation/common/OperationMatrix.hpp"

namespace sg
{

/**
 * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
 *
 * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data points
 *
 * This class defines an interface similar to OperationMultipleEval in order to support SIMD architectures
 * for datadriven task (multiple function evaluations, classification, regression). Target
 * architectures may be Intel SSE, Intel AVX, nVidia CUDA, OpenCL.
 */
class OperationMutlipleEvalVectorizedSP
{
protected:
	/// Pointer to the DataSet that should be evaluated on the grid
	DataMatrixSP* dataset_;
	/// Member to store the sparse grid's levels for better vectorization
	DataMatrixSP* level_;
	/// Member to store the sparse grid's indices for better vectorization
	DataMatrixSP* index_;

public:
	/**
	 * Constructor
	 */
	OperationMutlipleEvalVectorizedSP(DataMatrix& DataSet) {
		this->level_ = NULL;
		this->index_ = NULL;
	}

	/**
	 * Destructor
	 *
	 * cleans up level_ and index_ members
	 */
	virtual ~OperationMutlipleEvalVectorizedSP() {
		if (this->level_ != NULL)
			delete this->level_;

		if (this->index_ != NULL)
			delete this->index_;
	}

	/**
	 * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
	 *
	 * IMPORTANT REMARK:
	 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 	 * 	- data MUST a have even number of points AND it must be transposed
 	 *  - result MUST have the same size as data points that should be evaluated
	 *
	 * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void multVectorized(DataVectorSP& alpha, DataVectorSP& result) = 0;

	/**
	 * Multiplication of @f$B@f$ with vector @f$\alpha@f$
	 *
	 * IMPORTANT REMARK:
	 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 	 * 	- data MUST a have even number of points AND it must be transposed
 	 *  - result MUST have the same size as data points that should be evaluated
	 *
	 * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void multTransposeVectorized(DataVectorSP& source, DataVectorSP& result) = 0;

	/**
	 * rebuilds the DataMatrix for Level and Index in Derivatives
	 * needed for vectorization.
	 */
	virtual void rebuildLevelAndIndex() = 0;
};

}

#endif /* OPERATIONMULTIPLEEVALVECTORIZEDSP_HPP */
