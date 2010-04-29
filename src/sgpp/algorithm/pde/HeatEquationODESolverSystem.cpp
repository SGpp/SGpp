/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/HeatEquationODESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"

namespace sg
{

HeatEquationODESolverSystem::HeatEquationODESolverSystem(Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode)
{
	this->a = a;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;
	this->InnerGrid = NULL;
	this->alpha_inner = NULL;

	this->BoundaryUpdate = new DirichletUpdateVector(SparseGrid.getStorage());
	this->GridConverter = new DirichletGridConverter();

	this->OpLaplaceBound = SparseGrid.createOperationLaplace();
	this->OpMassBound = SparseGrid.createOperationLTwoDotProduct();

	// create the inner grid
	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

	//Create needed operations, on inner grid
	this->OpLaplaceInner = this->InnerGrid->createOperationLaplace();
	this->OpMassInner = this->InnerGrid->createOperationLTwoDotProduct();

	// right hand side if System
	this->rhs = new DataVector(1);
}

HeatEquationODESolverSystem::~HeatEquationODESolverSystem()
{
	delete this->OpLaplaceBound;
	delete this->OpMassBound;
	delete this->OpLaplaceInner;
	delete this->OpMassInner;

	delete this->BoundaryUpdate;
	delete this->GridConverter;
	if (this->InnerGrid != NULL)
	{
		delete this->InnerGrid;
	}
	if (this->alpha_inner != NULL)
	{
		delete this->alpha_inner;
	}
	delete this->rhs;
}

void HeatEquationODESolverSystem::mult(DataVector& alpha, DataVector& result)
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
	else
	{
		throw new algorithm_exception(" HeatEquationTimestepMatrix::mult : An unknown operation mode was specified!");
	}
}

DataVector* HeatEquationODESolverSystem::generateRHS()
{
	DataVector rhs_complete(this->alpha_complete->getSize());

	if (this->tOperationMode == "ExEul")
	{
		rhs_complete.setAll(0.0);

		DataVector temp(this->alpha_complete->getSize());

		applyMassMatrixComplete(*this->alpha_complete, temp);
		rhs_complete.add_parallel(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		rhs_complete.axpy_parallel(this->TimestepSize, temp);
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
		rhs_complete.add_parallel(temp);

		applyLOperatorComplete(*this->alpha_complete, temp);
		rhs_complete.axpy_parallel((0.5)*this->TimestepSize, temp);
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
		result_complete.add_parallel(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy_parallel((-1.0)*this->TimestepSize, temp);
	}
	else if (this->tOperationMode == "CrNic")
	{
		DataVector temp(alpha_bound.getSize());

		applyMassMatrixComplete(alpha_bound, temp);
		result_complete.add_parallel(temp);

		applyLOperatorComplete(alpha_bound, temp);
		result_complete.axpy_parallel((-0.5)*this->TimestepSize, temp);
	}
	else
	{
		throw new algorithm_exception("HeatEquationTimestepMatrix::generateRHS : An unknown operation mode was specified!");
	}
	rhs_complete.sub(result_complete);

	this->rhs->resize(this->alpha_inner->getSize());

	this->GridConverter->calcInnerCoefs(rhs_complete, *this->rhs);

	return this->rhs;
}

void HeatEquationODESolverSystem::applyMassMatrixComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());

	// Apply the mass matrix
	this->OpMassBound->mult(alpha, temp);

	result.add(temp);
}

void HeatEquationODESolverSystem::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());

	// Apply the laplace Operator rate
	this->OpLaplaceBound->mult(alpha, temp);
	result.axpy((-1.0)*this->a,temp);
}

void HeatEquationODESolverSystem::applyMassMatrixInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());

	// Apply the mass matrix
	this->OpMassInner->mult(alpha, temp);

	result.add(temp);
}

void HeatEquationODESolverSystem::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());

	// Apply the laplace Operator rate
	this->OpLaplaceInner->mult(alpha, temp);
	result.axpy((-1.0)*this->a,temp);
}

void HeatEquationODESolverSystem::finishTimestep()
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void HeatEquationODESolverSystem::startTimestep()
{
}

Grid* HeatEquationODESolverSystem::getGrid()
{
	return this->BoundGrid;
}

DataVector* HeatEquationODESolverSystem::getGridCoefficientsForCG()
{
	this->GridConverter->calcInnerCoefs(*this->alpha_complete, *this->alpha_inner);

	return this->alpha_inner;
}

DataVector* HeatEquationODESolverSystem::getGridCoefficients()
{
	return this->alpha_complete;
}

}
