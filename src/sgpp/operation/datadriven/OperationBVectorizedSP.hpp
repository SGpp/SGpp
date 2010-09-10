/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONBVECTORIZEDSP_HPP
#define OPERATIONBVECTORIZEDSP_HPP

#include "data/DataVectorSP.hpp"
#include "data/DataMatrixSP.hpp"

namespace sg
{

/**
 * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
 * 
 * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data points, then B is a @f$N\times m@f$ matrix, with
 * @f[ (B)_{i,j} = \varphi_i(x_j). @f]
 *
 * This class defines an interface similar to OperationB in order to support SIMD architectures
 * for datadriven task (multiple function evaluations, classification, regression). Target
 * architectures may be Intel SSE, Intel AVX, nVidia CUDA, OpenCL.
 *
 * This implementation uses single precision floating point numbers in order
 * to exploit consumer graphics boards.
 */
class OperationBVectorizedSP
{
public:
	/**
	 * Constructor
	 */
	OperationBVectorizedSP() {}

	/**
	 * Destructor
	 */
	virtual ~OperationBVectorizedSP() {}

	/**
	 * Multiplication of @f$B@f$ with vector @f$\alpha@f$
	 *
	 * IMPORTANT REMARK:
	 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 	 * 	- data MUST a have even number of points AND it must be transposed
 	 *  - result MUST have the same size as data points that should be evaluated
	 *
	 * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
	 * @param data vector, providing the data points x row-wise
	 * @param result the result vector of the matrix vector multiplication
	 *
	 * @return time needed to execute call
	 */
	virtual double multVectorized(DataVectorSP& alpha, DataMatrixSP& data, DataVectorSP& result) = 0;

	/**
	 * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
	 *
	 * IMPORTANT REMARK:
	 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 	 * 	- data MUST a have even number of points AND it must be transposed
 	 *  - result MUST have the same size as data points that should be evaluated
 	 *
	 * @param alpha vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
	 * @param data vector, providing the data points x row-wise
	 * @param result the result vector of the matrix vector multiplication
	 *
	 * @return time needed to execute call
	 */
	virtual double multTransposeVectorized(DataVectorSP& alpha, DataMatrixSP& data, DataVectorSP& result) = 0;

	/**
	 * rebuilds the DataMatrix for Level and Index in Derivatives
	 */
	virtual void rebuildLevelAndIndex() = 0;
};

}

#endif /* OPERATIONBVECTORIZEDSP_HPP */
