/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"
#include "exception/algorithm_exception.hpp"
using namespace sg::base;

namespace sg
{
namespace pde
{

PoissonEquationEllipticPDESolverSystemDirichlet::PoissonEquationEllipticPDESolverSystemDirichlet(Grid& SparseGrid, DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs)
{
	this->Laplace_Complete = this->BoundGrid->createOperationLaplace();
	this->Laplace_Inner = this->InnerGrid->createOperationLaplace();
}

PoissonEquationEllipticPDESolverSystemDirichlet::~PoissonEquationEllipticPDESolverSystemDirichlet()
{
	delete this->Laplace_Complete;
	delete this->Laplace_Inner;
}

void PoissonEquationEllipticPDESolverSystemDirichlet::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	Laplace_Inner->mult(alpha, result);
}

void PoissonEquationEllipticPDESolverSystemDirichlet::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	Laplace_Complete->mult(alpha, result);
}

}
}
