/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"
#include "exception/algorithm_exception.hpp"
#include "basis/operations_factory.hpp"

using namespace sg::base;
using namespace sg::GridOperationFactory;

namespace sg
{
namespace pde
{

PoissonEquationEllipticPDESolverSystemDirichlet::PoissonEquationEllipticPDESolverSystemDirichlet(Grid& SparseGrid, DataVector& rhs) : OperationEllipticPDESolverSystemDirichlet(SparseGrid, rhs)
{
	this->Laplace_Complete = createOperationLaplace(*this->BoundGrid);
	this->Laplace_Inner = createOperationLaplace(*this->InnerGrid);
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
