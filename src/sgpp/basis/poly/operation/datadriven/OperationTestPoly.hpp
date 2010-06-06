/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONTESTPOLY_HPP
#define OPERATIONTESTPOLY_HPP

#include "operation/common/OperationEval.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

/**
 * This class implements OperationTest for a grids with poly basis ansatzfunctions with
 *
 * @version $HEAD$
 */
class OperationTestPoly : public OperationTest
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the polynom's max. degree
	 */
	OperationTestPoly(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestPoly() {}

	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes);
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Poly Basis object
	SPolyBase base;
};

}

#endif /* OPERATIONTESTPOLY_HPP */
