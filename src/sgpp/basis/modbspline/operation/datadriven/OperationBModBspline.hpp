/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONBMODBSPLINE_HPP
#define OPERATIONBMODBSPLINE_HPP

#include "operation/datadriven/OperationB.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

/**
 * This class implements OperationB for a grid with modified Bspline basis functions
 *
 * @version $HEAD$
 */
class OperationBModBspline : public OperationB
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the Bspline's degree
	 */
	OperationBModBspline(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationBModBspline() {}

	virtual void mult(DataVector& alpha, DataVector& data, DataVector& result);
	virtual void multTranspose(DataVector& alpha, DataVector& data, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Mod Bspline Basis object
	SModBsplineBase base;
};

}

#endif /* OPERATIONBMODBSPLINE_HPP */
