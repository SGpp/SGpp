/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "operation/pde/OperationEllipticPDESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
using namespace sg::base;

namespace sg
{
namespace pde
{

OperationEllipticPDESolverSystem::OperationEllipticPDESolverSystem(Grid& SparseGrid, DataVector& rhs)
{
	this->BoundGrid = &SparseGrid;
	this->rhs = &rhs;

	this->numGridpointsComplete = this->BoundGrid->getSize();
}

OperationEllipticPDESolverSystem::~OperationEllipticPDESolverSystem()
{
}

size_t OperationEllipticPDESolverSystem::getNumGridPointsComplete()
{
	return this->numGridpointsComplete;
}

size_t OperationEllipticPDESolverSystem::getNumGridPointsInner()
{
	return this->numGridpointsInner;
}

}
}
