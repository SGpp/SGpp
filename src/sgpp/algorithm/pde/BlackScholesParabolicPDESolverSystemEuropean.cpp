/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/BlackScholesParabolicPDESolverSystemEuropean.hpp"
#include "exception/algorithm_exception.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "basis/operations_factory.hpp"
#include <cmath>

namespace sg
{
namespace finance
{

BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& mu,
			sg::base::DataVector& sigma, sg::base::DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
			bool bLogTransform, bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
			int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel)
{
	this->BoundGrid = &SparseGrid;
	this->alpha_complete = &alpha;

	this->alpha_complete_old = new sg::base::DataVector(*this->alpha_complete);
	this->alpha_complete_tmp = new sg::base::DataVector(*this->alpha_complete);

	this->InnerGrid = NULL;
	this->alpha_inner = NULL;
	this->tOperationMode = OperationMode;
	this->TimestepSize = TimestepSize;
	this->TimestepSize_old = TimestepSize;
	this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
	this->GridConverter = new sg::base::DirichletGridConverter();
	this->r = r;
	this->mus = &mu;
	this->sigmas = &sigma;
	this->rhos = &rho;
	this->BSalgoDims = this->BoundGrid->getAlgorithmicDimensions();
	this->nExecTimesteps = 1;

	// throw exception if grid dimensions not equal algorithmic dimensions
	if (this->BSalgoDims.size() != this->BoundGrid->getStorage()->dim())
	{
		throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : Number of algorithmic dimensions is not equal to the number of grid's dimensions.");
	}

	// test if number of dimensions in the coefficients match the numbers of grid dimensions (mu and sigma)
	if (this->BoundGrid->getStorage()->dim() != this->mus->getSize() || this->BoundGrid->getStorage()->dim() != this->sigmas->getSize())
	{
		throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : Dimension of mu and sigma parameters don't match the grid's dimensions!");
	}

	// test if number of dimensions in the coefficients match the numbers of grid dimensions (rho)
	if (this->BoundGrid->getStorage()->dim() != this->rhos->getNrows() || this->BoundGrid->getStorage()->dim() != this->rhos->getNcols())
	{
		throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : Row or col of rho parameter don't match the grid's dimensions!");
	}

	// test if all algorithmic dimensions are inside the grid's dimensions
	for (size_t i = 0; i < this->BSalgoDims.size(); i++)
	{
		if (this->BSalgoDims[i] >= this->BoundGrid->getStorage()->dim())
		{
			throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : Minimum one algorithmic dimension is not inside the grid's dimensions!");
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
			throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : There is minimum one doubled algorithmic dimension!");
		}
	}

	// build the coefficient vectors for the operations
	this->gammaCoef = new sg::base::DataMatrix(this->BSalgoDims.size(), this->BSalgoDims.size());
	this->deltaCoef = new sg::base::DataVector(this->BSalgoDims.size());

	// create the inner grid
	this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

	if (bLogTransform == false)
	{
		buildDeltaCoefficients();
		buildGammaCoefficients();

		//Create needed operations, on inner grid
		this->OpDeltaInner = sg::GridOperationFactory::createOperationDelta(*this->InnerGrid, *this->deltaCoef);
		this->OpGammaInner = sg::GridOperationFactory::createOperationGamma(*this->InnerGrid, *this->gammaCoef);
		// Create needed operations, on boundary grid
		this->OpDeltaBound = sg::GridOperationFactory::createOperationDelta(*this->BoundGrid, *this->deltaCoef);
		this->OpGammaBound = sg::GridOperationFactory::createOperationGamma(*this->BoundGrid, *this->gammaCoef);
	}
	// create needed operations that are different in case of a log-transformed Black-Scholoes equation
	else
	{
		buildDeltaCoefficientsLogTransform();
		buildGammaCoefficientsLogTransform();

		// operations on boundary grid
		this->OpDeltaBound = sg::GridOperationFactory::createOperationDeltaLog(*this->BoundGrid, *this->deltaCoef);
		this->OpGammaBound = sg::GridOperationFactory::createOperationGammaLog(*this->BoundGrid, *this->gammaCoef);
		//operations on inner grid
		this->OpDeltaInner = sg::GridOperationFactory::createOperationDeltaLog(*this->InnerGrid, *this->deltaCoef);
		this->OpGammaInner = sg::GridOperationFactory::createOperationGammaLog(*this->InnerGrid, *this->gammaCoef);
	}

	// Create operations, independent bLogTransform
	this->OpLTwoInner = sg::GridOperationFactory::createOperationLTwoDotProduct(*this->InnerGrid);
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

#ifdef HEDGE
	sg::base::BoundingBox* grid_bb = this->BoundGrid->getBoundingBox();
	sg::base::DimensionBoundary* myBoundaries = new sg::base::DimensionBoundary[grid_bb->getDimensions()];

	for (size_t d = 0; d < grid_bb->getDimensions(); d++)
	{
		double hedge_offset = (grid_bb->getIntervalWidth(d)-(grid_bb->getIntervalWidth(d)*HEDGE_WIDTH_PERCENT))/2.0;
		myBoundaries[d].leftBoundary = grid_bb->getBoundary(d).leftBoundary + hedge_offset;
		myBoundaries[d].rightBoundary = grid_bb->getBoundary(d).rightBoundary - hedge_offset;
		myBoundaries[d].bDirichletLeft = true;
		myBoundaries[d].bDirichletRight = true;
	}

	sg::base::BoundingBox* myHedgeBB = new sg::base::BoundingBox(grid_bb->getDimensions(), myBoundaries);
	// hedging
	myHedge = new sg::finance::Hedging(*myHedgeBB, 75, HEDGE_EPS);

	delete myHedgeBB;
	delete[] myBoundaries;
#endif
}

BlackScholesParabolicPDESolverSystemEuropean::~BlackScholesParabolicPDESolverSystemEuropean()
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
	delete this->alpha_complete_old;
	delete this->alpha_complete_tmp;

#ifdef HEDGE
	delete myHedge;
#endif
}

void BlackScholesParabolicPDESolverSystemEuropean::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

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

void BlackScholesParabolicPDESolverSystemEuropean::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

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

void BlackScholesParabolicPDESolverSystemEuropean::applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoBound->mult(alpha, temp);

	result.add(temp);
}

void BlackScholesParabolicPDESolverSystemEuropean::applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	sg::base::DataVector temp(alpha.getSize());

	result.setAll(0.0);

	// Apply the mass matrix
	this->OpLTwoInner->mult(alpha, temp);

	result.add(temp);
}

void BlackScholesParabolicPDESolverSystemEuropean::finishTimestep(bool isLastTimestep)
{
	// Replace the inner coefficients on the boundary grid
	this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

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
	this->numSumGridpointsInner += this->InnerGrid->getSize();
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

		// rebuild the inner grid + coefficients
		this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
	}

#ifdef HEDGE
	std::stringstream filename_ext;
	filename_ext << this->nExecTimesteps;
	this->nExecTimesteps++;

	myHedge->calc_hedging(*this->BoundGrid, *this->alpha_complete, filename_ext.str());
//	// Hedging (delta and gamma)
//	size_t nDims = this->BoundGrid->getStorage()->dim();
//	if (nDims <= 2)
//	{
//		std::stringstream filename_plot;
//		std::stringstream filename_delta;
//		std::stringstream filename_gamma;
//
//		filename_plot << "optionvalue_" << this->nExecTimesteps << ".plot";
//		filename_delta << "delta_" << this->nExecTimesteps << ".plot";
//		filename_gamma << "gamma_" << this->nExecTimesteps << ".plot";
//		this->nExecTimesteps++;
//
//		// get grid domain
//		base::BoundingBox* myBounds = this->BoundGrid->getBoundingBox();
//
//		// Plot option price
//		base::GridPrinter* myPrinter = new base::GridPrinter(*this->BoundGrid);
//		myPrinter->printGrid(*this->alpha_complete, filename_plot.str(), 100);
//		delete myPrinter;
//
//		// Plot delta
//		if (nDims == 1)
//		{
//			double left = myBounds->getBoundary(0).leftBoundary;
//			double right = myBounds->getBoundary(0).rightBoundary;
//			double diff_h = HEDGE_MESH_WIDTH;
//			std::ofstream deltaout;
//
//			base::OperationEval* myEval = GridOperationFactory::createOperationEval(*this->BoundGrid);
//			deltaout.open(filename_delta.str().c_str());
//
//			for (double pos = left + diff_h; pos < right - diff_h; pos+=(2.0*diff_h))
//			{
//				std::vector<double> point_left;
//				point_left.push_back(pos-diff_h);
//				std::vector<double> point_right;
//				point_right.push_back(pos+diff_h);
//
//				double res = (myEval->eval(*this->alpha_complete, point_right) - myEval->eval(*this->alpha_complete, point_left))/(2.0*diff_h);
//				deltaout << pos << " " << res << std::endl;
//			}
//
//			deltaout.close();
//			delete myEval;
//		}
//		else
//		{
//
//		}
//
//		// Plot gamma
//		if (nDims == 1)
//		{
//			double left = myBounds->getBoundary(0).leftBoundary;
//			double right = myBounds->getBoundary(0).rightBoundary;
//			double diff_h = HEDGE_MESH_WIDTH;
//			std::ofstream gammaout;
//
//			base::OperationEval* myEval = GridOperationFactory::createOperationEval(*this->BoundGrid);
//			gammaout.open(filename_gamma.str().c_str());
//
//			for (double pos = left + diff_h; pos < right - diff_h; pos+=(2.0*diff_h))
//			{
//				std::vector<double> point_left;
//				point_left.push_back(pos-diff_h);
//				std::vector<double> point_right;
//				point_right.push_back(pos+diff_h);
//				std::vector<double> point_middle;
//				point_middle.push_back(pos);
//
//				double res = (myEval->eval(*this->alpha_complete, point_right) - (2.0*myEval->eval(*this->alpha_complete, point_middle)) + myEval->eval(*this->alpha_complete, point_left))/(diff_h*diff_h);
//				gammaout << pos << " " << res << std::endl;
//			}
//
//			gammaout.close();
//			delete myEval;
//		}
//		else
//		{
//
//		}
//	}
#endif
}

void BlackScholesParabolicPDESolverSystemEuropean::startTimestep()
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

void BlackScholesParabolicPDESolverSystemEuropean::buildGammaCoefficients()
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

void BlackScholesParabolicPDESolverSystemEuropean::buildDeltaCoefficients()
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

void BlackScholesParabolicPDESolverSystemEuropean::buildGammaCoefficientsLogTransform()
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

void BlackScholesParabolicPDESolverSystemEuropean::buildDeltaCoefficientsLogTransform()
{
	size_t dim = this->BSalgoDims.size();

	for (size_t i = 0; i < dim; i++)
	{
		this->deltaCoef->set(i, this->mus->get(this->BSalgoDims[i])-(0.5*(this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[i]))));
	}
}

}
}
