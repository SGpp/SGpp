/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

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

	this->InnerGrid = NULL;
	this->alpha_inner = NULL;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->TimestepSize_old = TimestepSize;
	this->BoundaryUpdate = new DirichletUpdateVector(SparseGrid.getStorage());
	this->GridConverter = new DirichletGridConverter();
	this->a = a;
	this->theta = theta;
	this->sigma = sigma;

	// build the coefficient vectors for the operations
	//this->gammaCoef = new DataMatrix(SparseGrid.getStorage()->dim(), SparseGrid.getStorage()->dim());
	//this->deltaCoef = new DataVector(SparseGrid.getStorage()->dim());

	// create the inner grid
	//this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
	this->HWalgoDims = this->BoundGrid->getAlgorithmicDimensions();
	// create the inner grid
	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
	// Pass algorithmic dimensions to inner grid
	this->InnerGrid->setAlgorithmicDimensions(this->HWalgoDims);

	//	buildDeltaCoefficients();
	//	buildGammaCoefficients();

		//Create needed operations, on inner grid
		this->OpBInner = this->InnerGrid->createOperationLB();
		this->OpDInner = this->InnerGrid->createOperationLD();
		this->OpEInner = this->InnerGrid->createOperationLE();
	    this->OpFInner = this->InnerGrid->createOperationLF();
		// Create needed operations, on boundary grid
		this->OpBBound = this->BoundGrid->createOperationLB();
		this->OpDBound = this->BoundGrid->createOperationLD();
		this->OpEBound = this->BoundGrid->createOperationLE();
		this->OpFBound = this->BoundGrid->createOperationLF();

	// Create operations, independent bLogTransform
	this->OpLTwoInner = this->InnerGrid->createOperationLTwoDotProduct();
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
	delete this->OpBInner;
	delete this->OpDInner;
	delete this->OpEInner;
	delete this->OpFInner;
	delete this->OpLTwoInner;
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
	if (this->rhs != NULL)
	{
		delete this->rhs;
	}
	delete this->alpha_complete_old;
	delete this->alpha_complete_tmp;
}

void HullWhiteODESolverSystem::applyLOperatorComplete(DataVector& alpha, DataVector& result)
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

void HullWhiteODESolverSystem::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());
	result.setAll(0.0);
	    if (this->theta != 0.0)
			{
				this->OpBInner->mult(alpha, temp);
				result.axpy(1.0*this->theta, temp);
			}

		if (this->sigma != 0.0)
			{
				this->OpEInner->mult(alpha, temp);
				result.axpy((-1.0/2.0)*pow((this->sigma),2.0), temp);
			}

		if (this->a != 0.0)
			{
				this->OpFInner->mult(alpha, temp);
				result.axpy((-1.0)*this->a, temp);
			}


		this->OpDInner->mult(alpha, temp);
		result.sub(temp);

	//applyLOperatorComplete(alpha, result);
}

void HullWhiteODESolverSystem::applyMassMatrixComplete(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add(temp);
}

void HullWhiteODESolverSystem::applyMassMatrixInner(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoInner->mult(alpha, temp);

	result.add(temp);
	//applyMassMatrixComplete(alpha, result);
}

void HullWhiteODESolverSystem::finishTimestep(bool isLastTimestep)
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

	// Adjust the boundaries with the riskfree rate

		/*if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas")
		{
			this->BoundaryUpdate->multiplyBoundaryHullWhite(*this->alpha_complete,this->TimestepSize);
		}*/

	// add number of Gridpoints
	this->numSumGridpointsInner += this->InnerGrid->getSize();
	this->numSumGridpointsComplete += this->BoundGrid->getSize();

	if (this->useCoarsen == true && isLastTimestep == false)
	{
		size_t numCoarsen;

		// Coarsen the grid
		GridGenerator* myGeneratorCoarsen = this->BoundGrid->createGridGenerator();

		numCoarsen = myGeneratorCoarsen->getNumberOfRemoveablePoints();
		numCoarsen = static_cast<size_t>(((double)numCoarsen)*(this->coarsenPercent)/100.0);

		SurplusCoarseningFunctor* myCoarsenFunctor = new SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsenThreshold);

		for (size_t i = 0; i < this->numExecCoarsen; i++)
		{
			myGeneratorCoarsen->coarsen(myCoarsenFunctor, this->alpha_complete);
		}

		delete myGeneratorCoarsen;
		delete myCoarsenFunctor;

		// rebuild the inner grid + coefficients
		this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
	}
}

void HullWhiteODESolverSystem::startTimestep()
{/*
	   DataVector* factor = new DataVector(this->alpha_complete->getSize());
	// Adjust the boundaries with the riskfree rate
	   this->BoundaryUpdate->getfactor(*factor, this->TimestepSize);

		if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul")
		{
			this->BoundaryUpdate->multiplyBoundaryHullWhite(*this->alpha_complete,*factor);
		}*/

}
}
