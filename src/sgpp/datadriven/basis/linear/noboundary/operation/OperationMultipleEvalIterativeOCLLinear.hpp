/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVEOCLLINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVEOCLLINEAR_HPP

#include "datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "datadriven/operation/OCLKernels.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/tools/SGppStopwatch.hpp"

namespace sg
{
namespace parallel
{

/**
 * This class implements sg::base::OperationMultipleEval for a grids with linear basis ansatzfunctions without boundaries
 *
 * However in this case high efficient vector code (OpenCL) is generated
 * to implement a iterative OperationB version. In addition cache blocking is used
 * in order to assure a most efficient cache usage.
 *
 * IMPORTANT REMARK:
 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 * @li data MUST a have even number of points AND it must be transposed
 * @li result MUST have the same size as data points that should be evaluated
 */
class OperationMultipleEvalIterativeOCLLinear : public sg::base::OperationMultipleEvalVectorized
{
public:
	/**
	 * Construtor of OperationBLinear
	 *
	 * Within the construct sg::base::DataMatrix Level and sg::base::DataMatrix Index are set up.
	 * If the grid changes during your calculations and you don't want to create
	 * a new instance of this class, you have to call rebuildLevelAndIndex before
	 * doing any further mult or multTranspose calls.
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 * @param dataset dataset that should be evaluated
	 */
	OperationMultipleEvalIterativeOCLLinear(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset);

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalIterativeOCLLinear();

	virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result);

	virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result);

	virtual void rebuildLevelAndIndex();

protected:
	/// Pointer to the grid's gridstorage object
	sg::base::GridStorage* storage;
	/// sg::base::DataMatrix that contains a prepared Level matrix of all grid points (2^level)
	/// Timer object to handle time measurements
	sg::base::SGppStopwatch* myTimer;
	/// Object to access the OCL Kernel
	OCLKernels* myOCLKernels;
};

}
}

#endif /* OPERATIONMULTIPLEEVALITERATIVEOCLLINEAR_HPP */
