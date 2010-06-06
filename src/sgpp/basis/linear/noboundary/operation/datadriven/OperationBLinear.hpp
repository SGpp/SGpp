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

protected:
	/// Pointer to the grid's gridstorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONBLINEAR_HPP */
