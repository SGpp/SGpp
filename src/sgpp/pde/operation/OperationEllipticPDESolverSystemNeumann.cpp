/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/operation/OperationEllipticPDESolverSystemNeumann.hpp"
#include "base/exception/algorithm_exception.hpp"

namespace sg
{
namespace pde
{

OperationEllipticPDESolverSystemNeumann::OperationEllipticPDESolverSystemNeumann(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs) : OperationEllipticPDESolverSystem(SparseGrid, rhs)
{
}

OperationEllipticPDESolverSystemNeumann::~OperationEllipticPDESolverSystemNeumann()
{
}

void OperationEllipticPDESolverSystemNeumann::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	applyLOperator(alpha, result);
}

sg::base::DataVector* OperationEllipticPDESolverSystemNeumann::generateRHS()
{
	return this->rhs;
}

}
}

