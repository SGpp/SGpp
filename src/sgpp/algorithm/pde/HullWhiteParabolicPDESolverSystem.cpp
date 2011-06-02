/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi (qic@in.tum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "algorithm/pde/HullWhiteParabolicPDESolverSystem.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "grid/Grid.hpp"
#include "basis/operations_factory.hpp"

#include <cmath>

namespace sg
{
namespace finance
{

HullWhiteParabolicPDESolverSystem::HullWhiteParabolicPDESolverSystem(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, double sigma,
			double theta, double a, double TimestepSize, std::string OperationMode,
			bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
			int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel,
			int dim_HW)
{
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;

	this->alpha_complete_old = new sg::base::DataVector(this->alpha_complete->getSize());
	this->alpha_complete_old->setAll(0.0);
	this->alpha_complete_tmp = new sg::base::DataVector(this->alpha_complete->getSize());
	this->alpha_complete_tmp->setAll(0.0);
	this->alpha_complete_tmp->add(*this->alpha_complete);

	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->TimestepSize_old = TimestepSize;
	this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
	this->variableDiscountFactor = new VariableDiscountFactor(SparseGrid.getStorage(), dim_HW);
	this->a = a;
	this->theta = theta;
	this->sigma = sigma;
	this->HWalgoDims = this->BoundGrid->getAlgorithmicDimensions();

	// Create needed operations, on boundary grid
	this->OpBBound = sg::GridOperationFactory::createOperationLB(*this->BoundGrid);
	this->OpDBound = sg::GridOperationFactory::createOperationLD(*this->BoundGrid);
	this->OpEBound = sg::GridOperationFactory::createOperationLE(*this->BoundGrid);
	this->OpFBound = sg::GridOperationFactory::createOperationLF(*this->BoundGrid);

	// Create operations, independent bLogTransform
	this->OpLTwoBound = sg::GridOperationFactory::createOperationLTwoDotProduct(*this->BoundGrid);

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
	this->dim_r = dim_HW;
}

HullWhiteParabolicPDESolverSystem::~HullWhiteParabolicPDESolverSystem()
{
	delete this->OpBBound;
	delete this->OpDBound;
	delete this->OpEBound;
	delete this->OpFBound;
	delete this->OpLTwoBound;
	delete this->BoundaryUpdate;
	delete this->variableDiscountFactor;
	if (this->rhs != NULL)
	{
		delete this->rhs;
	}
	delete this->alpha_complete_old;
	delete this->alpha_complete_tmp;
}

void HullWhiteParabolicPDESolverSystem::applyLOperator(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

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

void HullWhiteParabolicPDESolverSystem::applyMassMatrix(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add(temp);
}

void HullWhiteParabolicPDESolverSystem::finishTimestep(bool isLastTimestep)
{
	sg::base::DataVector* factor = new sg::base::DataVector(this->alpha_complete->getSize());
	// Adjust the boundaries with the riskfree rate
	this->variableDiscountFactor->getDiscountFactor(*factor, this->TimestepSize);

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

void HullWhiteParabolicPDESolverSystem::startTimestep()
{
	   sg::base::DataVector* factor = new sg::base::DataVector(this->alpha_complete->getSize());
	// Adjust the boundaries with the riskfree rate
	   this->variableDiscountFactor->getDiscountFactor(*factor, this->TimestepSize);

		if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul")
		{
			this->BoundaryUpdate->multiplyBoundaryVector(*this->alpha_complete,*factor);
		}

}

}
}
