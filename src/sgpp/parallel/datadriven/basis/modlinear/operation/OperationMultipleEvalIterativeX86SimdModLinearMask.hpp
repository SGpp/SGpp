/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVEX86SIMDMODLINEARMASK_HPP
#define OPERATIONMULTIPLEEVALITERATIVEX86SIMDMODLINEARMASK_HPP

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"

namespace sg
{
namespace parallel
{

/**
 * This class implements sg::base::OperationMultipleEvalVectorized for a grids with modified
 * linear basis ansatzfunctions without boundaries
 *
 * However, in this case highly efficient vector code (AVX or SSE instructions) is generated
 * to implement an iterative OperationB version. In addition cache blocking is used
 * in order to assure a most efficient cache usage.
 *
 * IMPORTANT REMARK:
 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 * @li data MUST have a even number of points AND it must be transposed
 * @li result MUST have the same size as data points that should be evaluated
 */
class OperationMultipleEvalIterativeX86SimdModLinearMask : public sg::parallel::OperationMultipleEvalVectorized
{
public:
	/**
	 * Within the constructor, sg::base::DataMatrix Level and sg::base::DataMatrix Index are set up.
	 * If the grid changes during your calculations and you don't want to create
	 * a new instance of this class, you have to call rebuildLevelAndIndex before
	 * doing any further mult or multTranspose calls.
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 * @param dataset dataset that should be evaluated
	 */
	OperationMultipleEvalIterativeX86SimdModLinearMask(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset);

	virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result);

	virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result);

	virtual void rebuildLevelAndIndex();
};

}

}

#endif /* OPERATIONMULTIPLEEVALITERATIVEX86SIMDMODLINEARMASK_HPP */
