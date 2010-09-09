/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTMODPOLY_HPP
#define OPERATIONTESTMODPOLY_HPP

#include "operation/datadriven/OperationTest.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

/**
 * This class implements OperationTest for a grids with mod poly basis ansatzfunctions with
 *
 * @version $HEAD$
 */
class OperationTestModPoly : public OperationTest
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the polynom's max. degree
	 */
	OperationTestModPoly(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestModPoly() {}

	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes);
	virtual double testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues);
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Mod Poly Basis object
	SModPolyBase base;
};

}

#endif /* OPERATIONTESTMODPOLY_HPP */
