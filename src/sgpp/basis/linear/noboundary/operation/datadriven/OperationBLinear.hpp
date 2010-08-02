/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONBLINEAR_HPP
#define OPERATIONBLINEAR_HPP

#include "operation/datadriven/OperationB.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationB for a grids with linear basis ansatzfunctions without boundaries
 */
class OperationBLinear : public OperationB
{
public:
	/**
	 * Construtor of OperationBLinear
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationBLinear(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationBLinear() {}

	virtual void mult(DataVector& alpha, DataMatrix& data, DataVector& result);
	virtual void multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result);

	/**
	 * This method implements a highly data parallel and iterative evaluation of the sparse grid
	 * function at several datapoints provided in the data DataMatrix. Data parallelism is
	 * available when using the Intel compiler (icc). With the help of intrinsics high efficient
	 * SSE code is generated. You need a SSE3 compatible processor to execute this code.
	 *
	 * In order to use this routine you have to keep following points in mind:
	 * 	- data MUST a have even number of points AND it must be transposed
	 *  - result MUST have the same size as data points that should be evaluated
	 *
	 * @param Level DataMatrix that contains a prepared Level matrix of all grid points (2^level)
	 * @param Index DataMatrix that contains the indices off all grid points
	 * @param alpha DataVector containing the sparse grids ansatzfunctions' coefficients
	 * @param data DataMatrix containing the evaluation points; IMPROTANT: this Matrix must be transposed before calling this function
	 * @param result DataVector that will contain the result of the evaluation after execution
	 */
	void multTransposeIterative(DataMatrix& Level, DataMatrix& Index, DataVector& alpha, DataMatrix& data, DataVector& result);

	/**
	 * This method implements a highly data parallel and iterative transposed evaluation of the sparse grid
	 * function at several datapoints provided in the data DataMatrix. Data parallelism is
	 * available when using the Intel compiler (icc). With the help of intrinsics high efficient
	 * SSE code is generated. You need a SSE3 compatible processor to execute this code.
	 *
	 * In order to use this routine you have to keep following points in mind:
	 * 	- data MUST a have even number of points AND it must be transposed
	 *  - result MUST have the same size as data points that should be evaluated
	 *
	 * @param Level DataMatrix that contains a prepared Level matrix of all grid points (2^level)
	 * @param Index DataMatrix that contains the indices off all grid points
	 * @param source DataVector containing the sparse grids ansatzfunctions' coefficients
	 * @param data DataMatrix containing the evaluation points; IMPROTANT: this Matrix must be transposed before calling this function
	 * @param result DataVector that will contain the result of the evaluation after execution
	 */
	void multIterative(DataMatrix& Level, DataMatrix& Index, DataVector& source, DataMatrix& data, DataVector& result);

protected:
	/// Pointer to the grid's gridstorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONBLINEAR_HPP */
