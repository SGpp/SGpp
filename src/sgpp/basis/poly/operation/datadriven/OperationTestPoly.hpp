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
namespace datadriven
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
	 * @param storage the grid's sg::base::GridStorage object
	 * @param degree the polynom's max. degree
	 */
	OperationTestPoly(sg::base::GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestPoly() {}

	virtual double test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes);
	virtual double testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues);
	virtual double testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers);

protected:
	/// Pointer to sg::base::GridStorage object
	sg::base::GridStorage* storage;
	/// Poly Basis object
	SPolyBase base;
};

}
}

#endif /* OPERATIONTESTPOLY_HPP */
