/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/BlackScholesODESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include <cmath>

namespace sg
{

BlackScholesODESolverSystem::BlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu,
			DataVector& sigma, DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
			bool bLogTransform, bool useCoarsen, double coarsenThreshold, double coarsenPercent,
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
	this->r = r;
	this->mus = &mu;
	this->sigmas = &sigma;
	this->rhos = &rho;

	// build the coefficient vectors for the operations
	this->gammaCoef = new DataMatrix(SparseGrid.getStorage()->dim(), SparseGrid.getStorage()->dim());
	this->deltaCoef = new DataVector(SparseGrid.getStorage()->dim());

	// create the inner grid
	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

	if (bLogTransform == false)
	{
		buildDeltaCoefficients();
		buildGammaCoefficients();

		//Create needed operations, on inner grid
		this->OpDeltaInner = this->InnerGrid->createOperationDelta(*this->deltaCoef);
		this->OpGammaInner = this->InnerGrid->createOperationGamma(*this->gammaCoef);
		// Create needed operations, on boundary grid
		this->OpDeltaBound = this->BoundGrid->createOperationDelta(*this->deltaCoef);
		this->OpGammaBound = this->BoundGrid->createOperationGamma(*this->gammaCoef);
	}
	// create needed operations that are different in case of a log-transformed Black-Scholoes equation
	else
	{
		buildDeltaCoefficientsLogTransform();
		buildGammaCoefficientsLogTransform();

		// operations on boundary grid
		this->OpDeltaBound = this->BoundGrid->createOperationDeltaLog(*this->deltaCoef);
		this->OpGammaBound = this->BoundGrid->createOperationGammaLog(*this->gammaCoef);
		//operations on inner grid
		this->OpDeltaInner = this->InnerGrid->createOperationDeltaLog(*this->deltaCoef);
		this->OpGammaInner = this->InnerGrid->createOperationGammaLog(*this->gammaCoef);
	}

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

BlackScholesODESolverSystem::~BlackScholesODESolverSystem()
{
	delete this->OpDeltaBound;
	delete this->OpGammaBound;
	delete this->OpLTwoBound;
	delete this->OpDeltaInner;
	delete this->OpGammaInner;
	delete this->OpLTwoInner;
	delete this->gammaCoef;
	delete this->deltaCoef;
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
}

void BlackScholesODESolverSystem::applyLOperatorComplete(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpLTwoBound->mult(alpha, temp);
		result.axpy((-1.0)*this->r, temp);
	}

	// Apply the delta method
	this->OpDeltaBound->mult(alpha, temp);
	result.add(temp);

	// Apply the gamma method
	this->OpGammaBound->mult(alpha, temp);
	result.sub(temp);
}

void BlackScholesODESolverSystem::applyLOperatorInner(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the riskfree rate
	if (this->r != 0.0)
	{
		this->OpLTwoInner->mult(alpha, temp);
		result.axpy((-1.0)*this->r, temp);
	}

	// Apply the delta method
	this->OpDeltaInner->mult(alpha, temp);
	result.add(temp);

	// Apply the gamma method
	this->OpGammaInner->mult(alpha, temp);
	result.sub(temp);
}

void BlackScholesODESolverSystem::applyMassMatrixComplete(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add(temp);
}

void BlackScholesODESolverSystem::applyMassMatrixInner(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoInner->mult(alpha, temp);

	result.add(temp);
}

void BlackScholesODESolverSystem::finishTimestep(bool isLastTimestep)
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas")
		{
			this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}

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

		///////////////////////////////////////////////////
		// Start integrated refinement & coarsening
		///////////////////////////////////////////////////

//		size_t originalGridSize = this->BoundGrid->getStorage()->size();
//
//		// Coarsen the grid
//		GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();
//
//		size_t numRefines = myGenerator->getNumberOfRefinablePoints();
//		SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(this->alpha_complete, numRefines, this->coarsenThreshold);
//		//myGenerator->refineMaxLevel(myRefineFunc, this->numExecCoarsen);
//		myGenerator->refine(myRefineFunc);
//		this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
//		delete myRefineFunc;
//
//		size_t numCoarsen = myGenerator->getNumberOfRemoveablePoints();
//		SurplusCoarseningFunctor* myCoarsenFunctor = new SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, (this->coarsenThreshold/((double)this->numExecCoarsen)));
//		myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
//		delete myCoarsenFunctor;
//
//		delete myGenerator;

		///////////////////////////////////////////////////
		// End integrated refinement & coarsening
		///////////////////////////////////////////////////

		// rebuild the inner grid + coefficients
		this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
	}
}

void BlackScholesODESolverSystem::startTimestep()
{
	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul")
		{
			this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}
}

void BlackScholesODESolverSystem::buildGammaCoefficients()
{
	size_t dim = this->BoundGrid->getStorage()->dim();

	for (size_t i = 0; i < dim; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
			  this->gammaCoef->set(i, j, 0.5*((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get(i,j)));
			}
			else
			{
			  this->gammaCoef->set(i, j, ((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get(i,j)));
			}
		}
	}
}

void BlackScholesODESolverSystem::buildDeltaCoefficients()
{
	size_t dim = this->BoundGrid->getStorage()->dim();
	double covar_sum = 0.0;

	for (size_t i = 0; i < dim; i++)
	{
		covar_sum = 0.0;
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
				covar_sum += ((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get(i,j));
			}
			else
			{
				covar_sum += (0.5*((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get(i,j)));
			}
		}
		this->deltaCoef->set(i, this->mus->get(i)-covar_sum);
	}
}

void BlackScholesODESolverSystem::buildGammaCoefficientsLogTransform()
{
	size_t dim = this->BoundGrid->getStorage()->dim();

	for (size_t i = 0; i < dim; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
			  this->gammaCoef->set(i, j, 0.5*((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get(i,j)));
			}
			else
			{
			  this->gammaCoef->set(i, j, ((this->sigmas->get(i)*this->sigmas->get(j))*this->rhos->get(i,j)));
			}
		}
	}
}

void BlackScholesODESolverSystem::buildDeltaCoefficientsLogTransform()
{
	size_t dim = this->BoundGrid->getStorage()->dim();

	for (size_t i = 0; i < dim; i++)
	{
		this->deltaCoef->set(i, this->mus->get(i)-(0.5*(this->sigmas->get(i)*this->sigmas->get(i))));
	}
}

}
