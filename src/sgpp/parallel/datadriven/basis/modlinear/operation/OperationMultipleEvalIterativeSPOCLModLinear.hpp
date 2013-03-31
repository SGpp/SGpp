/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVESPOCLMODLINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVESPOCLMODLINEAR_HPP

#include "parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"
#include "parallel/datadriven/basis/common/OCLKernels.hpp"
#include "base/datatypes/DataVectorSP.hpp"

namespace sg
{

namespace parallel
{

/**
 * This class implements OperationMultipleEval for a grids with linear basis ansatzfunctions with modifications for boundaries
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
class OperationMultipleEvalIterativeSPOCLModLinear : public sg::parallel::OperationMultipleEvalVectorizedSP
{
public:
	/**
	 * Construtor of OperationBLinear
	 *
	 * Within the construct DataMatrix Level and DataMatrix Index are set up.
	 * If the grid changes during your calculations and you don't want to create
	 * a new instance of this class, you have to call rebuildLevelAndIndex before
	 * doing any further mult or multTranspose calls.
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 * @param dataset dataset that should be evaluated
	 */
	OperationMultipleEvalIterativeSPOCLModLinear(base::GridStorage* storage, base::DataMatrixSP* dataset);

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalIterativeSPOCLModLinear();

	virtual double multVectorized(base::DataVectorSP& alpha, base::DataVectorSP& result);

	virtual double multTransposeVectorized(base::DataVectorSP& source, base::DataVectorSP& result);

	virtual void rebuildLevelAndIndex();

protected:
	/// Object to access the OCL Kernel
	OCLKernels* myOCLKernels;
};

}

}

#endif /* OPERATIONMULTIPLEEVALITERATIVEOCLMODLINEAR_HPP */
