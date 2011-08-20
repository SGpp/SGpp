/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVEHYBRIDSSEOCLMODLINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVEHYBRIDSSEOCLMODLINEAR_HPP

#include "operation/datadriven/OperationMultipleEvalVectorized.hpp"
#include "basis/common/operation/datadriven/OCLKernels.hpp"
#include "grid/GridStorage.hpp"
#include "tools/common/SGppStopwatch.hpp"

namespace sg
{
namespace parallel
{

/**
 * This class implements sg::base::OperationMultipleEval for a grids with linear basis ansatzfunctions without boundaries
 * but with extrapolation towards the boundaries.
 *
 * However in this case high efficient vector code (OpenCL) is generated
 * to implement a iterative OperationB version. In addition cache blocking is used
 * in order to assure a most efficient cache usage. The CPU cores a are used by
 * an highly optimized SSE intrinsic implementation.
 *
 * IMPORTANT REMARK:
 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 * @li data MUST a have even number of points AND it must be transposed
 * @li result MUST have the same size as data points that should be evaluated
 */
class OperationMultipleEvalIterativeHybridSSEOCLModLinear : public sg::base::OperationMultipleEvalVectorized
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
	OperationMultipleEvalIterativeHybridSSEOCLModLinear(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset);

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalIterativeHybridSSEOCLModLinear();

	virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result);

	virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result);

	virtual void rebuildLevelAndIndex();

protected:
	/// Pointer to the grid's gridstorage object
	sg::base::GridStorage* storage;
	/// Timer object to handle time measurements
	sg::base::SGppStopwatch* myTimer;
	/// Object to access the OCL Kernel
	OCLKernels* myOCLKernels;
};

}

}

#endif /* OPERATIONMULTIPLEEVALITERATIVEHYBRIDSSEOCLMODLINEAR_HPP */
