/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALMODWAVELET_HPP
#define OPERATIONMULTIPLEEVALMODWAVELET_HPP

#include "operation/datadriven/OperationMultipleEval.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationMultipleEval for a grids with mod wavelet basis ansatzfunctions
 */
class OperationMultipleEvalModWavelet : public OperationMultipleEval
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationMultipleEvalModWavelet(GridStorage* storage, DataMatrix* dataset) : OperationMultipleEval(dataset) {
		this->storage = storage;
	}

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalModWavelet() {}

	virtual void mult(DataVector& alpha, DataVector& result);
	virtual void multTranspose(DataVector& source, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONMULTIPLEEVALMODWAVELET_HPP */
