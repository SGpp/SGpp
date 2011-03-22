/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVESPOCLLINEARSTRETCHED_HPP
#define OPERATIONMULTIPLEEVALITERATIVESPOCLLINEARSTRETCHED_HPP

#include "operation/datadriven/OperationMultipleEvalVectorizedSP.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OCLKernels.hpp"
#include "grid/GridStorage.hpp"
#include "tools/common/SGppStopwatch.hpp"

namespace sg
{

/**
 * This class implements OperationMultipleEvalSP for a grids with linearstretched basis ansatzfunctions without boundaries
 *
 * However in this case high efficient vector code (OpenCL) is generated
 * to implement a iterative OperationB version. In addition cache blocking is used
 * in order to assure a most efficient cache usage.
 *
 * IMPORTANT REMARK:
 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 * 	- data MUST a have even number of points AND it must be transposed
 *  - result MUST have the same size as data points that should be evaluated
 */
class OperationMultipleEvalIterativeSPOCLLinearStretched : public OperationMultipleEvalVectorizedSP
{
public:
	/**
	 * Construtor of OperationMultipleEvalLinearStretchedSP
	 *
	 * Within the construct DataMatrixSP Level and DataMatrixSP Index are set up.
	 * If the grid changes during your calculations and you don't want to create
	 * a new instance of this class, you have to call rebuildLevelAndIndex before
	 * doing any further mult or multTranspose calls.
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 * @param dataset dataset that should be evaluated
	 */
	OperationMultipleEvalIterativeSPOCLLinearStretched(GridStorage* storage, DataMatrixSP* dataset);

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalIterativeSPOCLLinearStretched();

	virtual double multVectorized(DataVectorSP& alpha,DataVectorSP& result);

	virtual double multTransposeVectorized(DataVectorSP& source, DataVectorSP& result);

	virtual void rebuildLevelAndIndex();

protected:
	/// Pointer to the grid's gridstorage object
	GridStorage* storage;
	/// Timer object to handle time measurements
	SGppStopwatch* myTimer;
	/// Object to access the OCL Kernel
	OCLKernels* myOCLKernels;
};

}

#endif /* OPERATIONMULTIPLEEVALITERATIVESPOCLLINEARSTRETCHED_HPP */
