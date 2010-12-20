/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi (qic@in.tum.de)

#include "algorithm/pde/HullWhiteODESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/Grid.hpp"
#include <cmath>

namespace sg
{

HullWhiteODESolverSystem::HullWhiteODESolverSystem(Grid& SparseGrid, DataVector& alpha, double sigma,
			double theta, double a, double TimestepSize, std::string OperationMode,
		    bool useCoarsen, double coarsenThreshold, double coarsenPercent,
			size_t numExecCoarsen)
{
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;

	this->alpha_complete_old = new DataVector(this->alpha_complete->getSize());
	this->alpha_complete_old->setAll(0.0);
	this->alpha_complete_tmp = new DataVector(this->alpha_complete->getSize());
	this->alpha_complete_tmp->setAll(0.0);
	this->alpha_complete_tmp->add(*this->alpha_complete);

	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->TimestepSize_old = TimestepSize;
	this->BoundaryUpdate = new DirichletUpdateVector(SparseGrid.getStorage());
	this->a = a;
	this->theta = theta;
	this->sigma = sigma;
	this->HWalgoDims = this->BoundGrid->getAlgorithmicDimensions();

	// Create needed operations, on boundary grid
	this->OpBBound = this->BoundGrid->createOperationLB();
	this->OpDBound = this->BoundGrid->createOperationLD();
	this->OpEBound = this->BoundGrid->createOperationLE();
	this->OpFBound = this->BoundGrid->createOperationLF();

	// Create operations, independent bLogTransform
	this->OpLTwoBound = this->BoundGrid->createOperationLTwoDotProduct();

	// right hand side if System
	this->rhs = NULL;

	// set coarsen settings
	this->useCoarsen = useCoarsen;
	this->coarsenThreshold = coarsenThreshold;
	this->coarsenPercent = coarsenPercent;
	this->numExecCoarsen = numExecCoarsen;

	// init Number of AverageGridPoins
	this->numSumGridpointsInner = 0;
	this->numSumGridpointsComplete = 0;
}

HullWhiteODESolverSystem::~HullWhiteODESolverSystem()
{
	delete this->OpBBound;
	delete this->OpDBound;
	delete this->OpEBound;
	delete this->OpFBound;
	delete this->OpLTwoBound;
	delete this->BoundaryUpdate;
	if (this->rhs != NULL)
	{
		delete this->rhs;
	}
	delete this->alpha_complete_old;
	delete this->alpha_complete_tmp;
}

void HullWhiteODESolverSystem::applyLOperator(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	if (this->theta != 0.0)
		{
			this->OpBBound->mult(alpha, temp);
			result.axpy(1.0*this->theta, temp);
		}

	if (this->sigma != 0.0)
		{
			this->OpEBound->mult(alpha, temp);
			result.axpy((-1.0/2.0)*pow((this->sigma),2.0), temp);
		}

	if (this->a != 0.0)
		{
			this->OpFBound->mult(alpha, temp);
			result.axpy((-1.0)*this->a, temp);
		}


	this->OpDBound->mult(alpha, temp);
	result.sub(temp);
}

void HullWhiteODESolverSystem::applyMassMatrix(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add(temp);
}

void HullWhiteODESolverSystem::finishTimestep(bool isLastTimestep)
{
	DataVector* factor = new DataVector(this->alpha_complete->getSize());
	// Adjust the boundaries with the riskfree rate
   this->BoundaryUpdate->getfactor(*factor, this->TimestepSize);

	if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas")
	{
		this->BoundaryUpdate->multiplyBoundaryVector(*this->alpha_complete,*factor);
	}

	// add number of Gridpoints
	this->numSumGridpointsInner += 0;
	this->numSumGridpointsComplete += this->BoundGrid->getSize();
}

void HullWhiteODESolverSystem::startTimestep()
{
	   DataVector* factor = new DataVector(this->alpha_complete->getSize());
	// Adjust the boundaries with the riskfree rate
	   this->BoundaryUpdate->getfactor(*factor, this->TimestepSize);

		if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul")
		{
			this->BoundaryUpdate->multiplyBoundaryVector(*this->alpha_complete,*factor);
		}

}

}
