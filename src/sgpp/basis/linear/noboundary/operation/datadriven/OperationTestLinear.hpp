/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTLINEAR_HPP
#define OPERATIONTESTLINEAR_HPP

#include "operation/datadriven/OperationTest.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationTest for a grids with linear basis ansatzfunctions without boundaries
 */
class OperationTestLinear : public OperationTest
{
public:
	/**
	 * Construtor of OperationTestLinear
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationTestLinear(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestLinear() {}

	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes);
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers);

protected:
	/// Pointer to the grid's gridstorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONTESTLINEAR_HPP */
