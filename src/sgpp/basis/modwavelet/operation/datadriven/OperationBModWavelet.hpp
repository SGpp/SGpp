/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONBMODWAVELET_HPP
#define OPERATIONBMODWAVELET_HPP

#include "operation/datadriven/OperationB.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationB for a grids with mod wavelet basis ansatzfunctions
 */
class OperationBModWavelet : public OperationB
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationBModWavelet(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationBModWavelet() {}

	virtual void mult(DataVector& alpha, DataMatrix& data, DataVector& result);
	virtual void multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONBMODWAVELET_HPP */
