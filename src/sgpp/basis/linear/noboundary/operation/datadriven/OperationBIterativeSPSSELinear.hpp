/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONBITERATIVESPSSELINEAR_HPP
#define OPERATIONBITERATIVESPSSELINEAR_HPP

#include "operation/datadriven/OperationBVectorizedSP.hpp"
#include "grid/GridStorage.hpp"
#include "tools/common/SGppStopwatch.hpp"

namespace sg
{

/**
 * This class implements OperationB for a grids with linear basis ansatzfunctions without boundaries
 *
 * However in this case high efficient vector code (SSE3 instructions) are generated
 * to implement a iterative OperationB version. In addition cache blocking is used
 * in order to assure a most efficient cache usage.
 *
 * IMPORTANT REMARK:
 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 * 	- data MUST a have even number of points AND it must be transposed
 *  - result MUST have the same size as data points that should be evaluated
 */
class OperationBIterativeSPSSELinear : public OperationBVectorizedSP
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
	 */
	OperationBIterativeSPSSELinear(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~OperationBIterativeSPSSELinear();

	virtual double multVectorized(DataVectorSP& alpha, DataMatrixSP& data, DataVectorSP& result);

	virtual double multTransposeVectorized(DataVectorSP& alpha, DataMatrixSP& data, DataVectorSP& result);

	virtual void rebuildLevelAndIndex();

protected:
	/// Pointer to the grid's gridstorage object
	GridStorage* storage;
	/// DataMatrix that contains a prepared Level matrix of all grid points (2^level)
	DataMatrixSP* Level;
	/// DataMatrix that contains the indices of all grid points
	DataMatrixSP* Index;
	/// Timer object to handle time measurements
	SGppStopwatch* myTimer;
};

}

#endif /* OPERATIONBITERATIVESPSSELINEAR_HPP */
