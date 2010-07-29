/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONBLINEARBOUNDARY_HPP
#define OPERATIONBLINEARBOUNDARY_HPP

#include "operation/datadriven/OperationB.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationB for a grids with linear basis ansatzfunctions
 * with boundaries
 *
 * @version $HEAD$
 */
class OperationBLinearBoundary : public OperationB
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GirdStorage object
	 */
	OperationBLinearBoundary(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationBLinearBoundary() {}

	virtual void mult(DataVector& alpha, DataMatrix& data, DataVector& result);
	virtual void multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONBLINEARBOUNDARY_HPP */
