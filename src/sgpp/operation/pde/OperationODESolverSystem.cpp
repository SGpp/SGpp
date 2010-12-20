/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "operation/pde/OperationODESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"

namespace sg
{

OperationODESolverSystem::OperationODESolverSystem()
{
	this->numSumGridpointsInner = 0;
	this->numSumGridpointsComplete = 0;
}

OperationODESolverSystem::~OperationODESolverSystem()
{
}

DataVector* OperationODESolverSystem::getGridCoefficients()
{
	return this->alpha_complete;
}

Grid* OperationODESolverSystem::getGrid()
{
	return this->BoundGrid;
}

void OperationODESolverSystem::setODESolver(std::string ode)
{
	this->tOperationMode = ode;
}

std::string OperationODESolverSystem::getODESolver()
{
	return this->tOperationMode;
}

void OperationODESolverSystem::setTimestepSize(double newTimestepSize)
{
	this->TimestepSize_old = this->TimestepSize;
	this->TimestepSize = newTimestepSize;
}

void OperationODESolverSystem::abortTimestep()
{
	*(this->alpha_complete) = *(this->alpha_complete_tmp);
}

void OperationODESolverSystem::saveAlpha()
{
	*(this->alpha_complete_old) = *(this->alpha_complete_tmp);
	*(this->alpha_complete_tmp) = *(this->alpha_complete);
}

size_t OperationODESolverSystem::getSumGridPointsComplete()
{
	return this->numSumGridpointsComplete;
}

size_t OperationODESolverSystem::getSumGridPointsInner()
{
	return this->numSumGridpointsInner;
}

}
