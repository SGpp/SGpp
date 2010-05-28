/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONEVALLINEAR_HPP
#define OPERATIONEVALLINEAR_HPP

#include "operation/common/OperationEval.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationEval for a grids with linear basis ansatzfunctions without boundaries
 */
class OperationEvalLinear : public OperationEval
{
public:
	/**
	 * Construtor of OperationEvalLinear
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationEvalLinear(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationEvalLinear() {}

	virtual double eval(DataVector& alpha, std::vector<double>& point);

protected:
	/// Pointer to the grid's gridstorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONEVALLINEAR_HPP */
