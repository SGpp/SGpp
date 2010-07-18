/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
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
}

OperationODESolverSystem::~OperationODESolverSystem()
{
}

void OperationODESolverSystem::mult(DataVector& alpha, DataVector& result)
{
	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixInner(alpha, result);
	}
	else if (this->tOperationMode == "ImEul")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrixInner(alpha, temp);
		result.add(temp);

		applyLOperatorInner(alpha, temp);
		result.axpy((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		result.setAll(0.0);

		DataVector temp(alpha.getSize());

		applyMassMatrixInner(alpha, temp);
		result.add(temp);

		applyLOperatorInner(alpha, temp);
		result.axpy((-0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "AdBas")
	{
		result.setAll(0.0);

		applyMassMatrixInner(alpha, result);
	}
	else
	{
		throw new algorithm_exception(" HeatEquationTimestepMatrix::mult : An unknown operation mode was specified!");
	}
}

DataVector* OperationODESolverSystem::generateRHS()
{
	DataVector rhs_complete(this->alpha_complete->getSize());

	if (this->tOperationMode == "ExEul")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		rhs_complete.axpy(this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "ImEul")
	{
		rhs_complete.setAll(0.0);

		applyMassMatrixComplete(*this->alpha_complete, rhs_complete);
	}
	else if (this->tOperationMode == "CrNic")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		rhs_complete.axpy((0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "AdBas")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);

		temp.mult((2.0)+this->TimestepSize/this->TimestepSize_old);

		DataVector temp_old(this->alpha_complete->getSize());
		applyMassMatrixComplete(*this->alpha_complete_old, temp_old);
		applyLOperatorComplete(*this->alpha_complete_old, temp_old);
		temp_old.mult(this->TimestepSize/this->TimestepSize_old);
		temp.sub(temp_old);

		rhs_complete.axpy((0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("HeatEquationTimestepMatrix::generateRHS : An unknown operation mode was specified!");
	}

	// Now we have the right hand side, lets apply the riskfree rate for the next timestep
	this->startTimestep();

	// Now apply the boundary ansatzfunctions to the inner ansatzfunctions
	DataVector result_complete(this->alpha_complete->getSize());
	DataVector alpha_bound(*this->alpha_complete);

	result_complete.setAll(0.0);

	this->BoundaryUpdate->setInnerPointsToZero(alpha_bound);

	// apply CG Matrix
	if (this->tOperationMode == "ExEul")
	{
		applyMassMatrixComplete(alpha_bound, result_complete);
	}
	else if (this->tOperationMode == "ImEul")
	{
		DataVector temp(alpha_bound.getSize());

		applyMassMatrixComplete(alpha_bound, temp);
		result_complete.add(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		DataVector temp(alpha_bound.getSize());

		applyMassMatrixComplete(alpha_bound, temp);
		result_complete.add(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy((-0.5)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "AdBas")
	{
		applyMassMatrixComplete(alpha_bound, result_complete);
	}
	else
	{
		throw new algorithm_exception("OperationODESolverSystem::generateRHS : An unknown operation mode was specified!");
	}
	rhs_complete.sub(result_complete);

	if (this->rhs != NULL)
	{
		delete this->rhs;
	}

	this->rhs = new DataVector(this->alpha_inner->getSize());
	this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

	return this->rhs;
}

void OperationODESolverSystem::finishTimestep()
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void OperationODESolverSystem::startTimestep()
{
}

DataVector* OperationODESolverSystem::getGridCoefficientsForCG()
{
	this->GridConverter->calcInnerCoefs(*this->alpha_complete, *this->alpha_inner);

	return this->alpha_inner;
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
	this->alpha_complete->setAll(0.0);
	this->alpha_complete->add(*this->alpha_complete_tmp);
}

void OperationODESolverSystem::saveAlpha() {
	this->alpha_complete_old->setAll(0.0);
	this->alpha_complete_old->add(*this->alpha_complete_tmp);
	this->alpha_complete_tmp->setAll(0.0);
	this->alpha_complete_tmp->add(*this->alpha_complete);
}

}
