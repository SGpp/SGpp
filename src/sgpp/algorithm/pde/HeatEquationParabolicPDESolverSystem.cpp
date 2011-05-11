/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/HeatEquationParabolicPDESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include "basis/operations_factory.hpp"
using namespace sg::base;
using namespace sg::GridOperationFactory;

namespace sg
{
namespace pde
{

HeatEquationParabolicPDESolverSystem::HeatEquationParabolicPDESolverSystem(Grid& SparseGrid, DataVector& alpha, double a, double TimestepSize, std::string OperationMode)
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

	this->OpLaplaceBound = createOperationLaplace(SparseGrid);
	this->OpMassBound = sg::GridOperationFactory::createOperationLTwoDotProduct(SparseGrid);

	// create the inner grid
	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

	//Create needed operations, on inner grid
	this->OpLaplaceInner = createOperationLaplace(*this->InnerGrid);
	this->OpMassInner = sg::GridOperationFactory::createOperationLTwoDotProduct(*this->InnerGrid);

	// right hand side if System
	this->rhs = new DataVector(1);
}

HeatEquationParabolicPDESolverSystem::~HeatEquationParabolicPDESolverSystem()
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

void HeatEquationParabolicPDESolverSystem::applyMassMatrixComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());

	// Apply the mass matrix
	this->OpMassBound->mult(alpha, temp);

	result.add(temp);
}

void HeatEquationParabolicPDESolverSystem::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());

	// Apply the laplace Operator rate
	this->OpLaplaceBound->mult(alpha, temp);
	result.axpy((-1.0)*this->a,temp);
}

void HeatEquationParabolicPDESolverSystem::applyMassMatrixInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());

	// Apply the mass matrix
	this->OpMassInner->mult(alpha, temp);

	result.add(temp);
}

void HeatEquationParabolicPDESolverSystem::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	DataVector temp(alpha.getSize());

	// Apply the laplace Operator rate
	this->OpLaplaceInner->mult(alpha, temp);
	result.axpy((-1.0)*this->a,temp);
}

void HeatEquationParabolicPDESolverSystem::finishTimestep(bool isLastTimestep)
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);
}

void HeatEquationParabolicPDESolverSystem::startTimestep()
{
}

}
}
