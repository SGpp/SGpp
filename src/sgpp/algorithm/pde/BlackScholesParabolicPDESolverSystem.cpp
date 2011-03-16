/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/BlackScholesParabolicPDESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include <cmath>
using namespace sg::base;

namespace sg
{

BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu,
			DataVector& sigma, DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
			bool bLogTransform, bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
			int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel)
{
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;

	this->alpha_complete_old = new DataVector(*this->alpha_complete);
	this->alpha_complete_tmp = new DataVector(*this->alpha_complete);

	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->TimestepSize_old = TimestepSize;
	this->BoundaryUpdate = new DirichletUpdateVector(SparseGrid.getStorage());
	this->r = r;
	this->mus = &mu;
	this->sigmas = &sigma;
	this->rhos = &rho;
	this->BSalgoDims = this->BoundGrid->getAlgorithmicDimensions();

	// throw exception if grid dimensions not equal algorithmic dimensions
	if (this->BSalgoDims.size() > this->BoundGrid->getStorage()->dim())
	{
		throw algorithm_exception("BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystemn : Number of algorithmic dimensions higher than the number of grid's dimensions.");
	}

	// test if number of dimensions in the coefficients match the numbers of grid dimensions (mu and sigma)
	if (this->BoundGrid->getStorage()->dim() != this->mus->getSize() || this->BoundGrid->getStorage()->dim() != this->sigmas->getSize())
	{
		throw algorithm_exception("BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem : Dimension of mu and sigma parameters don't match the grid's dimensions!");
	}

	// test if number of dimensions in the coefficients match the numbers of grid dimensions (rho)
	if (this->BoundGrid->getStorage()->dim() != this->rhos->getNrows() || this->BoundGrid->getStorage()->dim() != this->rhos->getNcols())
	{
		throw algorithm_exception("BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem : Row or col of rho parameter don't match the grid's dimensions!");
	}

	// test if all algorithmic dimensions are inside the grid's dimensions
	for (size_t i = 0; i < this->BSalgoDims.size(); i++)
	{
		if (this->BSalgoDims[i] >= this->BoundGrid->getStorage()->dim())
		{
			throw algorithm_exception("BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem : Minimum one algorithmic dimension is not inside the grid's dimensions!");
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
			throw algorithm_exception("BlackScholesParabolicPDESolverSystem::BlackScholesParabolicPDESolverSystem : There is minimum one doubled algorithmic dimension!");
		}
	}

	// build the coefficient vectors for the operations
	this->gammaCoef = new DataMatrix(this->BSalgoDims.size(), this->BSalgoDims.size());
	this->deltaCoef = new DataVector(this->BSalgoDims.size());

	if (bLogTransform == false)
	{
		buildDeltaCoefficients();
		buildGammaCoefficients();

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
	}

	this->OpLTwoBound = this->BoundGrid->createOperationLTwoDotProduct();

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

BlackScholesParabolicPDESolverSystem::~BlackScholesParabolicPDESolverSystem()
{
	delete this->OpDeltaBound;
	delete this->OpGammaBound;
	delete this->OpLTwoBound;
	delete this->gammaCoef;
	delete this->deltaCoef;
	delete this->BoundaryUpdate;
	if (this->rhs != NULL)
	{
		delete this->rhs;
	}
	delete this->alpha_complete_old;
	delete this->alpha_complete_tmp;
}

void BlackScholesParabolicPDESolverSystem::applyLOperator(DataVector& alpha, DataVector& result)
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

void BlackScholesParabolicPDESolverSystem::applyMassMatrix(DataVector& alpha, DataVector& result)
{
	DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add(temp);
}

void BlackScholesParabolicPDESolverSystem::finishTimestep(bool isLastTimestep)
{
#ifndef NOBOUNDARYDISCOUNT
	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas")
		{
			this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}
#endif

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
		GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

		//std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
		//std::cout << "Grid Size: " << originalGridSize << std::endl;

		if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine")
		{
			size_t numRefines = myGenerator->getNumberOfRefinablePoints();
			SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(this->alpha_complete, numRefines, this->refineThreshold);
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
			SurplusCoarseningFunctor* myCoarsenFunctor = new SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsenThreshold);
			myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
			delete myCoarsenFunctor;
		}

		delete myGenerator;

		///////////////////////////////////////////////////
		// End integrated refinement & coarsening
		///////////////////////////////////////////////////
	}
}

void BlackScholesParabolicPDESolverSystem::startTimestep()
{
#ifndef NOBOUNDARYDISCOUNT
	// Adjust the boundaries with the riskfree rate
	if (this->r != 0.0)
	{
		if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul")
		{
			this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0)*(this->r*this->TimestepSize))));
		}
	}
#endif
}

void BlackScholesParabolicPDESolverSystem::buildGammaCoefficients()
{
	size_t dim = this->BSalgoDims.size();

	for (size_t i = 0; i < dim; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
			  this->gammaCoef->set(i, j, 0.5*((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i],this->BSalgoDims[j])));
			}
			else
			{
			  this->gammaCoef->set(i, j, ((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i],this->BSalgoDims[j])));
			}
		}
	}
}

void BlackScholesParabolicPDESolverSystem::buildDeltaCoefficients()
{
	size_t dim = this->BSalgoDims.size();
	double covar_sum = 0.0;

	for (size_t i = 0; i < dim; i++)
	{
		covar_sum = 0.0;
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
				covar_sum += ((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i],this->BSalgoDims[j]));
			}
			else
			{
				covar_sum += (0.5*((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i],this->BSalgoDims[j])));
			}
		}
		this->deltaCoef->set(i, this->mus->get(this->BSalgoDims[i])-covar_sum);
	}
}

void BlackScholesParabolicPDESolverSystem::buildGammaCoefficientsLogTransform()
{
	size_t dim = this->BSalgoDims.size();

	for (size_t i = 0; i < dim; i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			// handle diagonal
			if (i == j)
			{
			  this->gammaCoef->set(i, j, 0.5*((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i],this->BSalgoDims[j])));
			}
			else
			{
			  this->gammaCoef->set(i, j, ((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i],this->BSalgoDims[j])));
			}
		}
	}
}

void BlackScholesParabolicPDESolverSystem::buildDeltaCoefficientsLogTransform()
{
	size_t dim = this->BSalgoDims.size();

	for (size_t i = 0; i < dim; i++)
	{
		this->deltaCoef->set(i, this->mus->get(this->BSalgoDims[i])-(0.5*(this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[i]))));
	}
}

}
