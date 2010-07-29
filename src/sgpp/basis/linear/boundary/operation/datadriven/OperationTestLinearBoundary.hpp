/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTLINEARBOUNDARY_HPP
#define OPERATIONTESTLINEARBOUNDARY_HPP

#include "operation/datadriven/OperationTest.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationTest for a grids with linear basis ansatzfunctions with
 * boundaries
 *
 * @version $HEAD$
 */
class OperationTestLinearBoundary : public OperationTest
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationTestLinearBoundary(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestLinearBoundary() {}

	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes);
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONTESTLINEARBOUNDARY_HPP */
