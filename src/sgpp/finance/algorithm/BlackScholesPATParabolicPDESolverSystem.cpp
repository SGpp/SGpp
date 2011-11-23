/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/algorithm/BlackScholesPATParabolicPDESolverSystem.hpp"
#include "base/exception/algorithm_exception.hpp"
#include "base/grid/generation/SurplusCoarseningFunctor.hpp"
#include "base/grid/generation/SurplusRefinementFunctor.hpp"
#include "pde/operation/PdeOpFactory.hpp"
#include <cmath>

namespace sg
{

namespace finance
{

BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& lambda,
			sg::base::DataMatrix& eigenvecs, sg::base::DataVector& mu_hat, double TimestepSize, std::string OperationMode,
			double dStrike, std::string option_type,
			bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
			int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel)
{
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;

	this->alpha_complete_old = new sg::base::DataVector(*this->alpha_complete);
	this->alpha_complete_tmp = new sg::base::DataVector(*this->alpha_complete);

	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->TimestepSize_old = TimestepSize;
	this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
	this->BSalgoDims = this->BoundGrid->getAlgorithmicDimensions();

	// set Eigenvalues, Eigenvector of covariance matrix and mu_hat
	this->lambda = new sg::base::DataVector(lambda);
	this->eigenvecs = new sg::base::DataMatrix(eigenvecs);
	this->mu_hat = new sg::base::DataVector(mu_hat);

	// throw exception if grid dimensions not equal algorithmic dimensions
	if (this->BSalgoDims.size() > this->BoundGrid->getStorage()->dim())
	{
		throw sg::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystemn : Number of algorithmic dimensions higher than the number of grid's dimensions.");
	}

	// test if number of dimensions in the coefficients match the numbers of grid dimensions (mu and sigma)
	if (this->BoundGrid->getStorage()->dim() != this->lambda->getSize())
	{
		throw sg::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem : Dimension of mu and sigma parameters don't match the grid's dimensions!");
	}

	// test if all algorithmic dimensions are inside the grid's dimensions
	for (size_t i = 0; i < this->BSalgoDims.size(); i++)
	{
		if (this->BSalgoDims[i] >= this->BoundGrid->getStorage()->dim())
		{
			throw sg::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem : Minimum one algorithmic dimension is not inside the grid's dimensions!");
		}
	}

	// test if there are double algorithmic dimensions
	std::vector<size_t> tempAlgoDims(this->BSalgoDims);
	for (size_t i = 0; i < this->BSalgoDims.size(); i++)
	{
		size_t dimCount = 0;
		for (size_t j = 0; j < tempAlgoDims.size(); j++)
		{
			if (this->BSalgoDims[i] == tempAlgoDims[j])
			{
				dimCount++;
			}
		}

		if (dimCount > 1)
		{
			throw sg::base::algorithm_exception("BlackScholesPATParabolicPDESolverSystem::BlackScholesPATParabolicPDESolverSystem : There is minimum one doubled algorithmic dimension!");
		}
	}

	// operations on boundary grid
	this->OpLaplaceBound = sg::op_factory::createOperationLaplace(*this->BoundGrid, *this->lambda);
	this->OpLTwoBound = sg::op_factory::createOperationLTwoDotProduct(*this->BoundGrid);

	// right hand side if System
	this->rhs = NULL;

	// set coarsen settings
	this->useCoarsen = useCoarsen;
	this->coarsenThreshold = coarsenThreshold;
	this->refineThreshold = refineThreshold;
	this->adaptSolveMode = adaptSolveMode;
	this->numCoarsenPoints = numCoarsenPoints;
	this->refineMode = refineMode;
	this->refineMaxLevel = refineMaxLevel;

	// init Number of AverageGridPoins
	this->numSumGridpointsInner = 0;
	this->numSumGridpointsComplete = 0;
}

BlackScholesPATParabolicPDESolverSystem::~BlackScholesPATParabolicPDESolverSystem()
{
	delete this->OpLaplaceBound;
	delete this->OpLTwoBound;
	delete this->BoundaryUpdate;
	if (this->rhs != NULL)
	{
		delete this->rhs;
	}
	delete this->alpha_complete_old;
	delete this->alpha_complete_tmp;
	delete this->lambda;
	delete this->eigenvecs;
	delete this->mu_hat;
}

void BlackScholesPATParabolicPDESolverSystem::applyLOperator(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());
	result.setAll(0.0);

	// Apply the Laplace operator
	this->OpLaplaceBound->mult(alpha, temp);
	result.axpy(-0.5, temp);
}

void BlackScholesPATParabolicPDESolverSystem::applyMassMatrix(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());
	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add(temp);
}

void BlackScholesPATParabolicPDESolverSystem::finishTimestep(bool isLastTimestep)
{
	// add number of Gridpoints
	this->numSumGridpointsInner += 0;
	this->numSumGridpointsComplete += this->BoundGrid->getSize();

	if (this->useCoarsen == true && isLastTimestep == false)
	{
		///////////////////////////////////////////////////
		// Start integrated refinement & coarsening
		///////////////////////////////////////////////////

		size_t originalGridSize = this->BoundGrid->getStorage()->size();

		// Coarsen the grid
		sg::base::GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

		//std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
		//std::cout << "Grid Size: " << originalGridSize << std::endl;

		if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine")
		{
			size_t numRefines = myGenerator->getNumberOfRefinablePoints();
			sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(this->alpha_complete, numRefines, this->refineThreshold);
			if (this->refineMode == "maxLevel")
			{
				myGenerator->refineMaxLevel(myRefineFunc, this->refineMaxLevel);
				this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
			}
			if (this->refineMode == "classic")
			{
				myGenerator->refine(myRefineFunc);
				this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
			}
			delete myRefineFunc;
		}

		if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine")
		{
			size_t numCoarsen = myGenerator->getNumberOfRemoveablePoints();
			sg::base::SurplusCoarseningFunctor* myCoarsenFunctor = new sg::base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsenThreshold);
			myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
			delete myCoarsenFunctor;
		}

		delete myGenerator;

		///////////////////////////////////////////////////
		// End integrated refinement & coarsening
		///////////////////////////////////////////////////
	}
}

void BlackScholesPATParabolicPDESolverSystem::startTimestep()
{
}

}

}
