/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVEHYBRIDX86SIMDOCLMODLINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVEHYBRIDX86SIMDOCLMODLINEAR_HPP

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "parallel/datadriven/basis/common/OCLKernels.hpp"
#include "parallel/tools/TwoPartitionAutoTuning.hpp"

namespace sg
{
namespace parallel
{

/**
 * This class implements sg::base::OperationMultipleEval for a grids with modified linear basis ansatzfunctions without boundaries
 * but with extrapolation towards the boundaries.
 *
 * However in this case highly efficient vector code (OpenCL and SSE/AVX intrinsics) is generated
 * to implement a iterative OperationB version. In addition cache blocking is used
 * in order to assure a most efficient cache usage.
 *
 * IMPORTANT REMARK:
 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 * @li data MUST a have even number of points AND it must be transposed
 * @li result MUST have the same size as data points that should be evaluated
 */
class OperationMultipleEvalIterativeHybridX86SimdOCLModLinear : public sg::parallel::OperationMultipleEvalVectorized
{
public:
	/**
	 * Construtor of sg::base::OperationMultipleEvalLinear
	 *
	 * Within the construct sg::base::DataMatrix Level and sg::base::DataMatrix Index are set up.
	 * If the grid changes during your calculations and you don't want to create
	 * a new instance of this class, you have to call rebuildLevelAndIndex before
	 * doing any further mult or multTranspose calls.
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 * @param dataset dataset that should be evaluated
	 */
	OperationMultipleEvalIterativeHybridX86SimdOCLModLinear(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset);

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalIterativeHybridX86SimdOCLModLinear();

	virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result);

	virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result);

	virtual void rebuildLevelAndIndex();

protected:
	/// Object to access the OCL Kernel
	OCLKernels* myOCLKernels;
	/// Autotuning object for mult routine
	sg::parallel::TwoPartitionAutoTuning* _tuningMult;
	/// Autotuning object for mult trans routine
	sg::parallel::TwoPartitionAutoTuning* _tuningMultTrans;
};

}

}

#endif /* OPERATIONMULTIPLEEVALITERATIVEHYBRIDX86SIMDOCLMODLINEAR_HPP */
