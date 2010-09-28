/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi (qic@in.tum.de)

#include "algorithm/pde/ModifiedBlackScholesODESolverSystem.hpp"
#include "algorithm/pde/BlackScholesODESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include <cmath>

namespace sg
{

ModifiedBlackScholesODESolverSystem::ModifiedBlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu,
			DataVector& sigma, DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
			bool bLogTransform, bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
			int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel)
: BlackScholesODESolverSystem(SparseGrid,
		alpha,
		mu,
		sigma,
		rho,
		r,
		TimestepSize,
		OperationMode,
		bLogTransform,
		useCoarsen,
		coarsenThreshold,
		adaptSolveMode,
		numCoarsenPoints,
		refineThreshold,
		refineMode,
		refineMaxLevel)
{
	this->OpFBound = this->BoundGrid->createOperationLF();
}

ModifiedBlackScholesODESolverSystem::~ModifiedBlackScholesODESolverSystem()
{
	delete this->OpFBound;
}

void ModifiedBlackScholesODESolverSystem::applyLOperatorComplete(DataVector& alpha, DataVector& result)
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

	//this->OpFBound->mult(alpha, temp);
	//result.axpy(0.06, temp);
	this->OpFBound->mult(alpha, temp);
	this->BoundaryUpdate->multiplyrBSHW(temp);
	//this->modifier->multiplyrBSHW(temp);
	result.add(temp);
}
void ModifiedBlackScholesODESolverSystem::finishTimestep(bool isLastTimestep)
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

void ModifiedBlackScholesODESolverSystem::startTimestep()
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
