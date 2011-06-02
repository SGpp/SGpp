/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "operation/pde/OperationParabolicPDESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include "basis/operations_factory.hpp"

namespace sg
{
namespace pde
{

OperationParabolicPDESolverSystem::OperationParabolicPDESolverSystem()
{
	this->numSumGridpointsInner = 0;
	this->numSumGridpointsComplete = 0;
}

OperationParabolicPDESolverSystem::~OperationParabolicPDESolverSystem()
{
}

sg::base::DataVector* OperationParabolicPDESolverSystem::getGridCoefficients()
{
	return this->alpha_complete;
}

sg::base::Grid* OperationParabolicPDESolverSystem::getGrid()
{
	return this->BoundGrid;
}

void OperationParabolicPDESolverSystem::setODESolver(std::string ode)
{
	this->tOperationMode = ode;
}

std::string OperationParabolicPDESolverSystem::getODESolver()
{
	return this->tOperationMode;
}

void OperationParabolicPDESolverSystem::setTimestepSize(double newTimestepSize)
{
	this->TimestepSize_old = this->TimestepSize;
	this->TimestepSize = newTimestepSize;
}

void OperationParabolicPDESolverSystem::abortTimestep()
{
	*(this->alpha_complete) = *(this->alpha_complete_tmp);
}

void OperationParabolicPDESolverSystem::saveAlpha()
{
	*(this->alpha_complete_old) = *(this->alpha_complete_tmp);
	*(this->alpha_complete_tmp) = *(this->alpha_complete);
}

size_t OperationParabolicPDESolverSystem::getSumGridPointsComplete()
{
	return this->numSumGridpointsComplete;
}

size_t OperationParabolicPDESolverSystem::getSumGridPointsInner()
{
	return this->numSumGridpointsInner;
}

void OperationParabolicPDESolverSystem::getGridCoefficientsForSC(sg::base::DataVector& Values)
{
	Values = *(this->alpha_complete);
	sg::base::OperationHierarchisation* myHierarchisation = sg::GridOperationFactory::createOperationHierarchisation(*BoundGrid);
	myHierarchisation->doDehierarchisation(Values);
	delete myHierarchisation;
}

}
}
