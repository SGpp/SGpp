/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONBMODPOLY_HPP
#define OPERATIONBMODPOLY_HPP

#include "operation/datadriven/OperationB.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

/**
 * This class implements OperationB for a grids with mod poly basis ansatzfunctions
 *
 * @version $HEAD$
 */
class OperationBModPoly : public OperationB
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the polynom's max. degree
	 */
	OperationBModPoly(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationBModPoly() {}

	virtual void mult(DataVector& alpha, DataMatrix& data, DataVector& result);
	virtual void multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Mod Poly Basis object
	SModPolyBase base;
};

}

#endif /* OPERATIONBMODPOLY_HPP */
