/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "operation/pde/OperationEllipticPDESolverSystemNeumann.hpp"
#include "exception/algorithm_exception.hpp"

namespace sg
{

OperationEllipticPDESolverSystemNeumann::OperationEllipticPDESolverSystemNeumann(Grid& SparseGrid, DataVector& rhs) : OperationEllipticPDESolverSystem(SparseGrid, rhs)
{
}

OperationEllipticPDESolverSystemNeumann::~OperationEllipticPDESolverSystemNeumann()
{
}

void OperationEllipticPDESolverSystemNeumann::mult(DataVector& alpha, DataVector& result)
{
	applyLOperator(alpha, result);
}

DataVector* OperationEllipticPDESolverSystemNeumann::generateRHS()
{
	return this->rhs;
}

}

