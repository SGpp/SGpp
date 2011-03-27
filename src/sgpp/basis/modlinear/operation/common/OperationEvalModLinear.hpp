/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONEVALMODLINEAR_HPP
#define OPERATIONEVELMODLINEAR_HPP

#include "operation/common/OperationEval.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationEval for a grids with mod linear basis ansatzfunctions with
 *
 * @version $HEAD$
 */
class OperationEvalModLinear : public OperationEval
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationEvalModLinear(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationEvalModLinear() {}

	virtual double eval(DataVector& alpha, std::vector<double>& point);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;

};

}

#endif /* OPERATIONEVELMODLINEAR_HPP */
