/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/algorithm/HestonParabolicPDESolverSystemEuroAmer.hpp"
#include "finance/application/HestonSolver.hpp"
#include "finance/application/BlackScholesSolver.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/ode/AdamsBashforth.hpp"
#include "solver/ode/VarTimestep.hpp"
#include "solver/ode/StepsizeControlH.hpp"
#include "solver/ode/StepsizeControlBDF.hpp"
#include "solver/ode/StepsizeControlEJ.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "base/grid/Grid.hpp"
#include "base/exception/application_exception.hpp"

#include "solver/sle/ConjugateGradients.hpp"

#include "base/operation/BaseOpFactory.hpp"
#include "pde/operation/PdeOpFactory.hpp"

#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <complex>
#include <limits>

using namespace sg::pde;
using namespace sg::solver;
using namespace sg::base;
using namespace std;

namespace sg
{
namespace finance
{

HestonSolver::HestonSolver(bool useLogTransform, bool usePAT) : ParabolicPDESolver()
{
	this->bStochasticDataAlloc = false;
	this->bGridConstructed = false;
	this->myScreen = NULL;
	this->useCoarsen = false;
	this->coarsenThreshold = 0.0;
	this->adaptSolveMode = "none";
	this->refineMode = "classic";
	this->numCoarsenPoints = -1;
	this->useLogTransform = useLogTransform;
	this->usePAT = usePAT;
	if (this->usePAT == true)
	{
		this->useLogTransform = true;
	}
	this->refineMaxLevel = 0;
	this->nNeededIterations = 0;
	this->dNeededTime = 0.0;
	this->staInnerGridSize = 0;
	this->finInnerGridSize = 0;
	this->avgInnerGridSize = 0;
	this->current_time = 0.0;
	this->tBoundaryType = "freeBoundaries";
}

HestonSolver::~HestonSolver()
{
	if (this->bStochasticDataAlloc)
	{
		delete this->thetas;
		delete this->kappas;
		delete this->volvols;
		delete this->hMatrix;
		//        delete this->eigval_covar;
		//        delete this->eigvec_covar;
		//        delete this->mu_hat;
	}
	if (this->myScreen != NULL)
	{
		delete this->myScreen;
	}
}

void HestonSolver::getGridNormalDistribution(DataVector& alpha, std::vector<double>& norm_mu, std::vector<double>& norm_sigma)
{
	if (this->bGridConstructed)
	{
		double value;
		StdNormalDistribution myNormDistr;
		double* s_coords = new double[this->dim];

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*(this->myBoundingBox));
			std::stringstream coordsStream(coords);

			for (size_t j = 0; j < this->dim; j++)
			{
				double tmp_load;
				coordsStream >> tmp_load;
				s_coords[j] = tmp_load;
			}

			value = 1.0;
			for (size_t j = 0; j < this->dim; j++)
			{
				if (this->useLogTransform == false)
				{
					value *= myNormDistr.getDensity(s_coords[j], norm_mu[j], norm_sigma[j]);
				}
				else
				{
					if (this->usePAT == true)
					{
						double inner_tmp = 0.0;
						for (size_t l = 0; l < dim; l++)
						{
							inner_tmp += this->eigvec_covar->get(j, l)*(s_coords[l]-(this->current_time*this->mu_hat->get(l)));
						}
						value *= myNormDistr.getDensity(exp(inner_tmp), norm_mu[j], norm_sigma[j]);
					}
					else
					{
						value *= myNormDistr.getDensity(exp(s_coords[j]), norm_mu[j], norm_sigma[j]);
					}
				}
			}

			alpha[i] = value;
		}
		delete[] s_coords;
	}
	else
	{
		throw new application_exception("HestonSolver::getGridNormalDistribution : The grid wasn't initialized before!");
	}
}

void HestonSolver::constructGrid(BoundingBox& BoundingBox, size_t level)
{
	this->dim = BoundingBox.getDimensions();

	if((dim%2) != 0)
		throw new application_exception("HestonSolver::constructGrid : The number of dimensions in the grid is not an even number! This doesn't correspond to an integer number of assets. The number of dimensions in the grid must be divisible by two.");

	this->numAssets = this->dim/2;
	this->levels = level;

	this->myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

	GridGenerator* myGenerator = this->myGrid->createGridGenerator();
	myGenerator->regular(this->levels);
	delete myGenerator;

	this->myBoundingBox = this->myGrid->getBoundingBox();
	this->myGridStorage = this->myGrid->getStorage();

	//std::string serGrid;
	//myGrid->serialize(serGrid);
	//std::cout << serGrid << std::endl;

	this->bGridConstructed = true;
}

void HestonSolver::refineInitialGridWithPayoff(DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance)
{
	size_t nRefinements = 0;

	this->dStrike = strike;
	this->payoffType = payoffType;

	if (this->useLogTransform == false)
	{
		if (this->bGridConstructed)
		{

			DataVector refineVector(alpha.getSize());

			if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put")
			{
				this->tBoundaryType = "Dirichlet";
				double tmp;
				double* dblFuncValues = new double[dim];
				double dDistance = 0.0;

				for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
				{
					std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*(this->myBoundingBox));
					std::stringstream coordsStream(coords);

					for (size_t j = 0; j < this->dim; j++)
					{
						coordsStream >> tmp;

						dblFuncValues[j] = tmp;
					}

					tmp = 0.0;
					for (size_t j = 0; j < this->dim; j++)
					{
						tmp += dblFuncValues[j];
					}

					if (payoffType == "std_euro_call")
					{
						dDistance = fabs(((tmp/static_cast<double>(this->dim))-strike));
					}
					if (payoffType == "std_euro_put" || payoffType == "std_amer_put")
					{
						dDistance = fabs((strike-(tmp/static_cast<double>(this->dim))));
					}

					if (dDistance <= dStrikeDistance)
					{
						refineVector[i] = dDistance;
						nRefinements++;
					}
					else
					{
						refineVector[i] = 0.0;
					}
				}

				delete[] dblFuncValues;

				SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&refineVector, nRefinements, 0.0);

				this->myGrid->createGridGenerator()->refine(myRefineFunc);

				delete myRefineFunc;

				alpha.resize(this->myGridStorage->size());

				// reinit the grid with the payoff function
				initGridWithPayoff(alpha, strike, payoffType);
			}
			else
			{
				throw new application_exception("HestonSolver::refineInitialGridWithPayoff : An unsupported payoffType was specified!");
			}
		}
		else
		{
			throw new application_exception("HestonSolver::refineInitialGridWithPayoff : The grid wasn't initialized before!");
		}
	}
}

void HestonSolver::refineInitialGridWithPayoffToMaxLevel(DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance, size_t maxLevel)
{
	size_t nRefinements = 0;

	this->dStrike = strike;
	this->payoffType = payoffType;

	if (this->useLogTransform == false)
	{
		if (this->bGridConstructed)
		{

			DataVector refineVector(alpha.getSize());

			if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put")
			{
				double tmp;
				double* dblFuncValues = new double[dim];
				double dDistance = 0.0;

				this->tBoundaryType = "Dirichlet";

				for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
				{
					std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
					std::stringstream coordsStream(coords);

					for (size_t j = 0; j < this->dim; j++)
					{
						coordsStream >> tmp;

						dblFuncValues[j] = tmp;
					}

					tmp = 0.0;
					for (size_t j = 0; j < this->dim; j++)
					{
						tmp += dblFuncValues[j];
					}

					if (payoffType == "std_euro_call")
					{
						dDistance = fabs(((tmp/static_cast<double>(this->dim))-strike));
					}
					if (payoffType == "std_euro_put" || payoffType == "std_amer_put")
					{
						dDistance = fabs((strike-(tmp/static_cast<double>(this->dim))));
					}

					if (dDistance <= dStrikeDistance)
					{
						refineVector[i] = dDistance;
						nRefinements++;
					}
					else
					{
						refineVector[i] = 0.0;
					}
				}

				delete[] dblFuncValues;

				SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&refineVector, nRefinements, 0.0);

				this->myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

				delete myRefineFunc;

				alpha.resize(this->myGridStorage->size());

				// reinit the grid with the payoff function
				initGridWithPayoff(alpha, strike, payoffType);
			}
			else
			{
				throw new application_exception("HestonSolver::refineInitialGridWithPayoffToMaxLevel : An unsupported payoffType was specified!");
			}
		}
		else
		{
			throw new application_exception("HestonSolver::refineInitialGridWithPayoffToMaxLevel : The grid wasn't initialized before!");
		}
	}
}

void HestonSolver::setStochasticData(DataVector& thetas_arg, DataVector& kappas_arg, DataVector& volvols_arg, DataMatrix& hMatrix_arg, double r)
{
	this->thetas = new sg::base::DataVector(thetas_arg);
	this->kappas = new sg::base::DataVector(kappas_arg);
	this->volvols = new sg::base::DataVector(volvols_arg);
	this->hMatrix = new sg::base::DataMatrix(hMatrix_arg);
	this->r = r;

	// calculate eigenvalues, eigenvectors and mu_hat from stochastic data for PAT
	// (not required yet)
	//    size_t mydim = this->thetas->getSize();
	//    this->eigval_covar = new sg::base::DataVector(mydim);
	//    this->eigvec_covar = new sg::base::DataMatrix(mydim,mydim);
	//    this->mu_hat = new sg::base::DataVector(mydim);

	//PAT stuff (not required yet)
	//    for (size_t i = 0; i < mydim; i++)
	//      {
	//        double tmp = 0.0;
	//        for (size_t j = 0; j < mydim; j++)
	//          {
	//            tmp += ((this->mus->get(j) - (0.5*this->sigmas->get(j)*this->sigmas->get(j))) * this->eigvec_covar->get(j, i));
	//          }
	//        this->mu_hat->set(i, tmp);
	//      }

	bStochasticDataAlloc = true;
}

void HestonSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	//    if (this->bGridConstructed && this->bStochasticDataAlloc)
	//      {
	//        Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
	//        SLESolver* myCG = NULL;
	//        OperationParabolicPDESolverSystem* myHestonSystem = NULL;
	//
	//        if (this->tBoundaryType == "Dirichlet")
	//          {
	//           myCG = new BiCGStab(maxCGIterations, epsilonCG);
	//           myHestonSystem = new HestonParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->thetas, *this->volvols, *this->kappas, *this->rhos, this->r, timestepsize, "ExEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//          }
	//        else
	//          {
	//
	//          }
	//
	//        SGppStopwatch* myStopwatch = new SGppStopwatch();
	//        this->staInnerGridSize = getNumberInnerGridPoints();
	//
	//        std::cout << "Using Explicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
	//        myStopwatch->start();
	//        myEuler->solve(*myCG, *myHestonSystem, true, verbose);
	//        this->dNeededTime = myStopwatch->stop();
	//
	//        std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
	//        std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;
	//
	//        std::cout << "Average Grid size: " << static_cast<double>(myHestonSystem->getSumGridPointsComplete())/static_cast<double>(numTimesteps) << std::endl;
	//        std::cout << "Average Grid size (Inner): " << static_cast<double>(myHestonSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;
	//
	//        if (this->myScreen != NULL)
	//          {
	//            std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
	//            this->myScreen->writeEmptyLines(2);
	//          }
	//
	//        this->finInnerGridSize = getNumberInnerGridPoints();
	//        this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myHestonSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps))+0.5);
	//        this->nNeededIterations = myEuler->getNumberIterations();
	//
	//        delete myHestonSystem;
	//        delete myCG;
	//        delete myEuler;
	//        delete myStopwatch;
	//
	//        this->current_time += (static_cast<double>(numTimesteps)*timestepsize);
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::solveExplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
	//      }
}

void HestonSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	//    if (this->bGridConstructed && this->bStochasticDataAlloc)
	//      {
	//        Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
	//        SLESolver* myCG = NULL;
	//        OperationParabolicPDESolverSystem* myBSSystem = NULL;
	//
	//        if (this->tBoundaryType == "Dirichlet")
	//          {
	//            if (this->usePAT == true)
	//              {
	//                myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
	//#ifdef _OPENMP
	//                myBSSystem = new BlackScholesPATParabolicPDESolverSystemEuroAmerParallelOMP(*this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat, timestepsize, "ImEul", this->dStrike, this->payoffType, this->r, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//#else
	//                myBSSystem = new BlackScholesPATParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat, timestepsize, "ImEul", this->dStrike, this->payoffType, this->r, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//#endif
	//              }
	//            else
	//              {
	//                myCG = new BiCGStab(maxCGIterations, epsilonCG);
	//#ifdef _OPENMP
	//                myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//#else
	//                myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//#endif
	//              }
	//          }
	//        else
	//          {
	//            if (this->usePAT == true)
	//              {
	//                myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
	//                myBSSystem = new BlackScholesPATParabolicPDESolverSystem(*this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//              }
	//            else
	//              {
	//                myCG = new BiCGStab(maxCGIterations, epsilonCG);
	//                myBSSystem = new BlackScholesParabolicPDESolverSystem(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//              }
	//          }
	//
	//        SGppStopwatch* myStopwatch = new SGppStopwatch();
	//        this->staInnerGridSize = getNumberInnerGridPoints();
	//
	//        std::cout << "Using Implicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
	//        myStopwatch->start();
	//        myEuler->solve(*myCG, *myBSSystem, true, verbose);
	//        this->dNeededTime = myStopwatch->stop();
	//
	//        std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
	//        std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;
	//
	//        std::cout << "Average Grid size: " << static_cast<double>(myBSSystem->getSumGridPointsComplete())/static_cast<double>(numTimesteps) << std::endl;
	//        std::cout << "Average Grid size (Inner): " << static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;
	//
	//        if (this->myScreen != NULL)
	//          {
	//            std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
	//            this->myScreen->writeEmptyLines(2);
	//          }
	//
	//        this->finInnerGridSize = getNumberInnerGridPoints();
	//        this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps))+0.5);
	//        this->nNeededIterations = myEuler->getNumberIterations();
	//
	//        delete myBSSystem;
	//        delete myCG;
	//        delete myEuler;
	//        delete myStopwatch;
	//
	//        this->current_time += (static_cast<double>(numTimesteps)*timestepsize);
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::solveImplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
	//      }
}

void HestonSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul)
{
	if (this->bGridConstructed && this->bStochasticDataAlloc)
	{
		SLESolver* myCG = NULL;
		OperationParabolicPDESolverSystem* myHestonSystem = NULL;

		if (this->tBoundaryType == "Dirichlet")
		{
			myCG = new BiCGStab(maxCGIterations, epsilonCG);

			//			GridStorage* myGridStorage = this->myGrid->getStorage();
			//
			//			// determine the number of grid points for both grids
			//			size_t numTotalGridPoints = myGridStorage->size();
			//			size_t numInnerGridPoints = myGridStorage->getNumInnerPoints();
			//
			//			size_t numNonZeroInner = 0;
			//			for (size_t i = 0; i < numTotalGridPoints; i++)
			//			{
			//				GridIndex* curPoint = (*myGridStorage)[i];
			//				if(curPoint->isInnerPoint() && alpha.get(i) != 0)
			//					numNonZeroInner++;
			//			}

			myHestonSystem = new HestonParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->thetas, *this->volvols, *this->kappas, *this->hMatrix, this->r, timestepsize, "CrNic", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
		}
		else
		{
		}

		SGppStopwatch* myStopwatch = new SGppStopwatch();
		this->staInnerGridSize = getNumberInnerGridPoints();

		size_t numCNSteps;
		size_t numIESteps;

		numCNSteps = numTimesteps;
		if (numTimesteps > NumImEul)
		{
			numCNSteps = numTimesteps - NumImEul;
		}
		numIESteps = NumImEul;

		Euler* myEuler = new Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
		CrankNicolson* myCN = new CrankNicolson(numCNSteps, timestepsize, this->myScreen);

		myStopwatch->start();
		if (numIESteps > 0)
		{
			std::cout << "Using Implicit Euler to solve " << numIESteps << " timesteps:" << std::endl;
			myHestonSystem->setODESolver("ImEul");
			myEuler->solve(*myCG, *myHestonSystem, false, false);
		}
		myHestonSystem->setODESolver("CrNic");
		std::cout << "Using Crank Nicolson to solve " << numCNSteps << " timesteps:" << std::endl << std::endl << std::endl << std::endl;
		myCN->solve(*myCG, *myHestonSystem, true, false);
		this->dNeededTime = myStopwatch->stop();

		std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
		std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

		std::cout << "Average Grid size: " << static_cast<double>(myHestonSystem->getSumGridPointsComplete())/static_cast<double>(numTimesteps) << std::endl;
		std::cout << "Average Grid size (Inner): " << static_cast<double>(myHestonSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		this->finInnerGridSize = getNumberInnerGridPoints();
		this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myHestonSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps))+0.5);
		this->nNeededIterations = myEuler->getNumberIterations() + myCN->getNumberIterations();

		delete myHestonSystem;
		delete myCG;
		delete myCN;
		delete myEuler;
		delete myStopwatch;

		this->current_time += (static_cast<double>(numTimesteps)*timestepsize);
	}
	else
	{
		throw new application_exception("HestonSolver::solveCrankNicolson : A grid wasn't constructed before or stochastic parameters weren't set!");
	}
}


void HestonSolver::solveAdamsBashforth(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose)
{
	//    ODESolver* myODESolver = new AdamsBashforth(numTimesteps, timestepsize, myScreen);
	//    HestonSolver::solveX(numTimesteps, timestepsize, maxCGIterations, epsilonCG, alpha, verbose, myODESolver, "AdBas");
	//    delete myODESolver;
}

void HestonSolver::solveSC(std::string Solver, size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose)
{
	std::string tmp;
	float epsilon = 0.001;
	float sc = 1;
	float gamma = 0.5;
	ODESolver* myODESolver;
	std::istringstream iss(Solver);
	if(Solver[2] == '2') {
		getline(iss,tmp,':');
		getline(iss,tmp,':');
		std::istringstream qwe(tmp);
		qwe >> epsilon;
		iss >> gamma;
		std::cout<<"2 " << "AdBas"<<", "  <<"CrNic"  << " Epsilon: "<< epsilon << " Gamma: "  << gamma   << std::endl;
		myODESolver = new VarTimestep("AdBas","CrNic",numTimesteps, timestepsize, epsilon, myScreen, gamma);

	} else if (Solver[2] == 'H') {
		getline(iss,tmp,':');
		getline(iss,tmp,':');
		std::istringstream qwe(tmp);
		qwe >> epsilon;
		iss >> gamma;
		std::cout<< "H "  <<"CrNic"  << " Epsilon: "<< epsilon << " Gamma: "  << gamma   << std::endl;
		myODESolver = new StepsizeControlH("CrNic",numTimesteps, timestepsize, epsilon, myScreen, gamma);


	}else if (Solver[2] == 'I') {
		getline(iss,tmp,':');
		getline(iss,tmp,':');
		std::istringstream qwe(tmp);
		qwe >> epsilon;
		getline(iss,tmp,':');
		std::istringstream qwe2(tmp);
		qwe >> sc;
		iss >> gamma;
		std::cout << "I "   << " Epsilon: "<< epsilon<< " SC: "<<sc<< " Gamma: "  << gamma   << std::endl;
		myODESolver = new StepsizeControlEJ("CrNic",numTimesteps, timestepsize, epsilon, sc,  myScreen, gamma);

	} else std::cerr << "HestonSolver::solveSC(): Unknown Stepsize Control #" << Solver[3] << "#" << Solver << std::endl;

	HestonSolver::solveX(numTimesteps, timestepsize, maxCGIterations, epsilonCG, alpha, verbose, myODESolver, "CrNic");
	delete myODESolver;
}

void HestonSolver::solveSCAC(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose)
{
	ODESolver* myODESolver = new VarTimestep("AdBasC","CrNic",numTimesteps, timestepsize, epsilon, myScreen, -1);
	HestonSolver::solveX(numTimesteps, timestepsize, maxCGIterations, epsilonCG, alpha, verbose, myODESolver, "CrNic");
	delete myODESolver;
}

void HestonSolver::solveSCH(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose)
{
	ODESolver* myODESolver = new StepsizeControlH("CrNic",numTimesteps, timestepsize, epsilon, myScreen, 0.9);
	HestonSolver::solveX(numTimesteps, timestepsize, maxCGIterations, epsilonCG, alpha, verbose, myODESolver, "CrNic");
	delete myODESolver;
}

void HestonSolver::solveSCBDF(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose)
{
	ODESolver* myODESolver = new StepsizeControlBDF(numTimesteps, timestepsize, epsilon, myScreen);
	HestonSolver::solveX(numTimesteps, timestepsize, maxCGIterations, epsilonCG, alpha, verbose, myODESolver, "SCBDF");
	delete myODESolver;
}

void HestonSolver::solveSCEJ(size_t numTimesteps, double timestepsize, double epsilon, double myAlpha, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose)
{
	ODESolver* myODESolver = new StepsizeControlEJ("CrNic",numTimesteps, timestepsize, epsilon, myAlpha,  myScreen, 0.5);
	HestonSolver::solveX(numTimesteps, timestepsize, maxCGIterations, epsilonCG, alpha, verbose, myODESolver, "SCEJ");
	delete myODESolver;
}

void HestonSolver::solveX(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, void *myODESolverV, std::string Solver)
{
	//    ODESolver *myODESolver = (ODESolver *)myODESolverV;
	//    if (this->bGridConstructed && this->bStochasticDataAlloc)
	//      {		BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
	//        OperationParabolicPDESolverSystem* myBSSystem = NULL;
	//
	//        if (this->tBoundaryType == "Dirichlet")
	//          {
	//#ifdef _OPENMP
	//            myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, Solver, this->dStrike, this->payoffType, this->useLogTransform, false, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//#else
	//            myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, Solver, this->dStrike, this->payoffType, this->useLogTransform, false, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//#endif
	//          }
	//        else
	//          {
	//            myBSSystem = new BlackScholesParabolicPDESolverSystem(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, Solver, this->dStrike, this->payoffType, this->useLogTransform, false, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//          }
	//
	//        SGppStopwatch* myStopwatch = new SGppStopwatch();
	//        this->staInnerGridSize = getNumberInnerGridPoints();
	//
	//        myStopwatch->start();
	//        myODESolver->solve(*myCG, *myBSSystem, false, verbose);
	//        this->dNeededTime = myStopwatch->stop();
	//
	//        if (this->myScreen != NULL)
	//          {
	//            std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
	//            this->myScreen->writeEmptyLines(2);
	//          }
	//
	//        this->finInnerGridSize = getNumberInnerGridPoints();
	//        this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps))+0.5);
	//        this->nNeededIterations = myODESolver->getNumberIterations();
	//
	//        delete myBSSystem;
	//        delete myCG;
	//
	//        this->current_time += (static_cast<double>(numTimesteps)*timestepsize);
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::solveX : A grid wasn't constructed before or stochastic parameters weren't set!");
	//      }
}


void HestonSolver::initGridWithPayoff(DataVector& alpha, double strike, std::string payoffType)
{
	this->dStrike = strike;
	this->payoffType = payoffType;

	if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put")
	{
		this->tBoundaryType = "Dirichlet";
	}

	if (this->useLogTransform)
	{
		if (this->usePAT)
		{
			initPATTransformedGridWithPayoff(alpha, strike, payoffType);
		}
		else
		{
			initLogTransformedGridWithPayoff(alpha, strike, payoffType);
		}
	}
	else
	{
		initCartesianGridWithPayoff(alpha, strike, payoffType);

		//		sg::base::GridStorage* myGridStorage = myGridStorage;

		// determine the number of grid points for both grids
		size_t numTotalGridPoints = myGridStorage->size();
		size_t numInnerGridPoints = myGridStorage->getNumInnerPoints();

		size_t numNonZeroInner = 0;
		for (size_t i = 0; i < numTotalGridPoints; i++)
		{
			sg::base::GridIndex* curPoint = (*myGridStorage)[i];
			if(curPoint->isInnerPoint() && alpha.get(i) != 0)
				numNonZeroInner++;
		}

		int k=0;
		k++;

	}
}

double HestonSolver::get1DEuroCallPayoffValue(double assetValue, double strike)
{
	if (assetValue <= strike)
	{
		return 0.0;
	}
	else
	{
		return assetValue - strike;
	}
}

double HestonSolver::getAnalyticSolution1D(double stock, bool isCall, double t, double vola, double r, double strike)
{
	StdNormalDistribution myStdNDis;

	double dOne = (log((stock/strike)) + ((r + (vola*vola*0.5))*(t)))/(vola*sqrt(t));
	double dTwo = dOne - (vola*sqrt(t));

	if (isCall)
	{
		return (stock*myStdNDis.getCumulativeDensity(dOne)) - (strike*myStdNDis.getCumulativeDensity(dTwo)*(exp((-1.0)*r*t)));
	}
	else
	{
		return (strike*myStdNDis.getCumulativeDensity(dTwo*(-1.0))*(exp((-1.0)*r*t))) - (stock*myStdNDis.getCumulativeDensity(dOne*(-1.0)));
	}
}


void HestonSolver::solve1DAnalytic(std::vector< std::pair<double, double> >& premiums, double minStock, double maxStock, double StockInc, double strike, double t, bool isCall)
{
	//    if (bStochasticDataAlloc)
	//      {
	//        double stock = 0.0;
	//        double vola = this->sigmas->get(0);
	//
	//        for (stock = minStock; stock <= maxStock; stock += StockInc)
	//          {
	//            double prem = getAnalyticSolution1D(stock, isCall, t, vola, this->r, strike);
	//            premiums.push_back(std::make_pair(stock, prem));
	//          }
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::solve1DAnalytic : Stochastic parameters weren't set!");
	//      }
}

void HestonSolver::print1DAnalytic(std::vector< std::pair<double, double> >& premiums, std::string tfilename)
{
	typedef std::vector< std::pair<double, double> > printVector;
	std::ofstream fileout;

	fileout.open(tfilename.c_str());
	for(printVector::iterator iter = premiums.begin(); iter != premiums.end(); iter++)
	{
		fileout << iter->first << " " << iter->second << " " << std::endl;
	}
	fileout.close();
}


void HestonSolver::getAnalyticAlpha1D(DataVector& alpha_analytic, double strike, double t, std::string payoffType, bool hierarchized)
{
	//    double coord;
	//
	//    if(dim!=1)
	//      {
	//        throw new application_exception("HestonSolver::getAnalyticAlpha1D : A grid wasn't constructed before!");
	//      }
	//    if (!this->bGridConstructed)
	//      {
	//        throw new application_exception("HestonSolver::getAnalyticAlpha1D : function only available for dim = 1!");
	//      }
	//
	//    // compute values of analytic solution on given grid
	//    for (size_t i = 0; i < this->myGridStorage->size(); i++)
	//      {
	//        std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
	//        std::stringstream coordsStream(coords);
	//        coordsStream >> coord;
	//        if(useLogTransform)
	//          {
	//            coord = exp(coord);
	//          }
	//        if (payoffType == "std_euro_call")
	//          {
	//            alpha_analytic[i] = getAnalyticSolution1D(coord, true, t, this->sigmas->get(0), this->r, strike);
	//          }
	//        else if (payoffType == "std_euro_put")
	//          {
	//            alpha_analytic[i] = getAnalyticSolution1D(coord, false, t, this->sigmas->get(0), this->r, strike);
	//          }
	//      }
	//
	//    if(hierarchized)
	//      {
	//        // hierarchize computed values
	//        OperationHierarchisation* myHier = sg::op_factory::createOperationHierarchisation(*this->myGrid);
	//        myHier->doHierarchisation(alpha_analytic);
	//
	//        delete myHier;
	//      }
}


void HestonSolver::evaluate1DAnalyticCuboid(sg::base::DataVector& AnalyticOptionPrices, sg::base::DataMatrix& EvaluationPoints, double strike, double vola, double r, double t, bool isCall)
{
	size_t n = EvaluationPoints.getNrows();

	if (AnalyticOptionPrices.getSize() != n)
	{
		throw new sg::base::application_exception("PDESolver::evaluate1DAnalyticCuboid : The size of the price vector doesn't match the size of the evaluation points' vector!");
	}

	for(size_t k=0; k<n; k++)
	{
		double x = EvaluationPoints.get(k,0); // get first coordinate

		if(this->useLogTransform)
		{
			x = exp(x);
		}
		double price = getAnalyticSolution1D(x, isCall, t, vola, r, strike);
		AnalyticOptionPrices.set(k,price);
	}
}

std::vector<size_t> HestonSolver::getAlgorithmicDimensions()
{
	return this->myGrid->getAlgorithmicDimensions();
}

void HestonSolver::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims)
{
	if (this->tBoundaryType == "freeBoundaries")
	{
		this->myGrid->setAlgorithmicDimensions(newAlgoDims);
	}
	else
	{
		throw new application_exception("HestonSolver::setAlgorithmicDimensions : Set algorithmic dimensions is only supported when choosing option type all!");
	}
}

void HestonSolver::initScreen()
{
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Heston Solver, 1.0.0", "TUM (C) 2009-2010, by Sam Maurus");
	this->myScreen->writeStartSolve("Multidimensional Heston Solver");
}

void HestonSolver::setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, size_t refineMaxLevel, int numCoarsenPoints, double coarsenThreshold, double refineThreshold)
{
	this->useCoarsen = true;
	this->coarsenThreshold = coarsenThreshold;
	this->refineThreshold = refineThreshold;
	this->refineMaxLevel = refineMaxLevel;
	this->adaptSolveMode = adaptSolveMode;
	this->refineMode = refineMode;
	this->numCoarsenPoints = numCoarsenPoints;
}

void HestonSolver::printPayoffInterpolationError2D(DataVector& alpha, std::string tFilename, size_t numTestpoints, double strike)
{
	if (this->useLogTransform == false)
	{
		if (this->bGridConstructed)
		{
			if (this->myGrid->getStorage()->getBoundingBox()->getDimensions() == 2)
			{
				if (numTestpoints < 2)
					numTestpoints = 2;

				double dInc = (2.0*strike)/static_cast<double>(numTestpoints-1);

				double dX = 0.0;
				double dY = 2*strike;

				std::ofstream file;
				file.open(tFilename.c_str());

				OperationEval* myEval = sg::op_factory::createOperationEval(*this->myGrid);

				for (size_t i = 0; i < numTestpoints; i++)
				{
					std::vector<double> point;

					point.push_back(dX);
					point.push_back(dY);

					double result = myEval->eval(alpha, point);

					file << std::scientific << std::setprecision( 16 ) << dX << " " << dY << " " << result << std::endl;

					dX += dInc;
					dY -= dInc;
				}

				delete myEval;

				file.close();
			}
		}
		else
		{
			throw new application_exception("HestonSolver::getPayoffInterpolationError : A grid wasn't constructed before!");
		}
	}
}

size_t HestonSolver::getGridPointsAtMoney(std::string payoffType, double strike, double eps)
{
	size_t nPoints = 0;

	if (this->useLogTransform == false)
	{
		if (this->bGridConstructed)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				bool isAtMoney = true;
				DataVector coords(this->dim);
				this->myGridStorage->get(i)->getCoordsBB(coords, *this->myBoundingBox);

				// The stock prices come from the coordinates of the 0th, 2nd, 4th, 6th... dimensions. The other dimensions are the volatilities...they don't count for adding up stock prices.
				// So...in order to determine whether the point is at the money or not, we have to add up only the 0th, 2nd, 4th, 6th...dimension values.
				if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put")
				{
					double stockPriceSum = 0.0;
					for(size_t j=0;j<this->dim;j=j+2)
					{
						stockPriceSum += coords[j];
					}

					if ( ((stockPriceSum/static_cast<double>(this->dim)) < (strike-eps)) || ((stockPriceSum/static_cast<double>(this->dim)) > (strike+eps)) )
					{
						isAtMoney = false;
					}
				}
				else
				{
					throw new application_exception("HestonSolver::getGridPointsAtMoney : An unknown payoff-type was specified!");
				}

				if (isAtMoney == true)
				{
					nPoints++;
				}
			}
		}
		else
		{
			throw new application_exception("HestonSolver::getGridPointsAtMoney : A grid wasn't constructed before!");
		}
	}

	return nPoints;
}

void HestonSolver::initCartesianGridWithPayoff(DataVector& alpha, double strike, std::string payoffType)
{
	double tmp;

	BlackScholesSolver* bsSolver = new BlackScholesSolver();

	if (this->bGridConstructed)
	{
		std::ofstream fileout;

		//		fileout.open("innerpoints.txt");

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);

			GridIndex* curPoint = (*myGridStorage)[i];

			std::stringstream coordsStream(coords);
			double* dblFuncValues = new double[dim];

			for (size_t j = 0; j < this->dim; j++)
			{
				coordsStream >> tmp;

				dblFuncValues[j] = tmp;
			}

			// print the values to a file

			if (payoffType == "std_euro_call")
			{
				// So now we have in alpha[i] the standard payoff value for the point S_i, v_i.
				// What we now want to do is solve the exact BS equation for every point S_i, v_i;
				//				alpha[i] = bsSolver->getAnalyticSolution1D(dblFuncValues[0], true, 0.05, pow(dblFuncValues[1], 2.0), this->r, this->dStrike);

				tmp = 0.0;
				for (size_t j = 0; j < numAssets; j++)
				{
					tmp += dblFuncValues[2*j];
				}

				//				tmp = 0.0;
				//				for (size_t j = 0; j < dim; j++)
				//				{
				//					tmp += dblFuncValues[j];
				//				}

				if(!curPoint->isInnerPoint() && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).leftBoundary)
				{
					// Use the Black-Scholes solution for this boundary
					alpha[i] = bsSolver->getAnalyticSolution1D(dblFuncValues[0], true, 1.0, pow(dblFuncValues[1], 2.0), this->r, this->dStrike);
					//					alpha[i] = (this->myBoundingBox->getBoundary(0).rightBoundary-strike)/(this->myBoundingBox->getBoundary(0).rightBoundary)*dblFuncValues[0];
				}
				else if(!curPoint->isInnerPoint() && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).rightBoundary)
				{
					// Dirichlet condition when v -> inf is that U = S
					alpha[i] = dblFuncValues[0]; //(this->myBoundingBox->getBoundary(0).rightBoundary-strike)/(this->myBoundingBox->getBoundary(0).rightBoundary)*dblFuncValues[0];
				}
				else if(!curPoint->isInnerPoint() && dblFuncValues[0] == this->myBoundingBox->getBoundary(0).rightBoundary)
				{
					// Set boundary to be the linear function at s_max
					double normalPayoff = std::max<double>(((tmp/static_cast<double>(numAssets))-strike), 0.0);
					double vRange = this->myBoundingBox->getBoundary(1).rightBoundary - this->myBoundingBox->getBoundary(1).leftBoundary;
					double sPayoffDiff = dblFuncValues[0] - normalPayoff;
					alpha[i] = normalPayoff + ((dblFuncValues[1] - this->myBoundingBox->getBoundary(1).leftBoundary) /vRange)*sPayoffDiff;
				}
				else
					alpha[i] = std::max<double>(((tmp/static_cast<double>(numAssets))-strike), 0.0);

				//				fileout << i << " " << dblFuncValues[0] << " " << dblFuncValues[1] << " " << alpha[i] << " inner point " << curPoint->isInnerPoint() << "" << std::endl;
			}
			else if (payoffType == "std_euro_put" || payoffType == "std_amer_put")
			{
				tmp = 0.0;
				for (size_t j = 0; j < dim; j++)
				{
					tmp += dblFuncValues[j];
				}
				alpha[i] = std::max<double>(strike-((tmp/static_cast<double>(dim))), 0.0);
			}
			else
			{
				throw new application_exception("HestonSolver::initCartesianGridWithPayoff : An unknown payoff-type was specified!");
			}

			delete[] dblFuncValues;
		}

		//		fileout.close();

		// determine the number of grid points for both grids
		size_t numTotalGridPoints = myGridStorage->size();
		size_t numInnerGridPoints = myGridStorage->getNumInnerPoints();

		size_t numNonZeroInner = 0;
		for (size_t i = 0; i < numTotalGridPoints; i++)
		{
			sg::base::GridIndex* curPoint = (*myGridStorage)[i];
			if(curPoint->isInnerPoint() && alpha.get(i) != 0)
				numNonZeroInner++;
		}

		int k=0;
		k++;

		OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("HestonSolver::initCartesianGridWithPayoff : A grid wasn't constructed before!");
	}
}

void HestonSolver::initLogTransformedGridWithPayoff(DataVector& alpha, double strike, std::string payoffType)
{
	double tmp;

	BlackScholesSolver* bsSolver = new BlackScholesSolver();

	if (this->bGridConstructed)
	{
		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
			std::stringstream coordsStream(coords);
			double* dblFuncValues = new double[dim];

			GridIndex* curPoint = (*myGridStorage)[i];

			for (size_t j = 0; j < this->dim; j++)
			{
				coordsStream >> tmp;

				dblFuncValues[j] = tmp;
			}

			if (payoffType == "std_euro_call")
			{
				// Here we have to be a little careful in comparison to the Black-Scholes solver.
				// In the Black-Scholes solver, every dimension represents an asset price, so determining the payoff value at this particular point
				// is simply a case of adding up all the dimension (asset) values and inserting that into the payoff function.
				// In the Heston model, however, only every second dimension represents an asset price. Every other dimension represents a variance, which
				// we MUST NOT use to determine the payoff value.
				// Thus, in oder to determine the payoff value, we have to only iterate over the first, third, fifth... dimensions.

				tmp = 0.0;
				for (size_t j = 0; j < numAssets; j++)
				{
					tmp += exp(dblFuncValues[2*j]);
				}

				if(!curPoint->isInnerPoint() && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).leftBoundary)
				{
					// Use the Black-Scholes solution for this boundary
					alpha[i] = bsSolver->getAnalyticSolution1D(tmp, true, 1.0, pow(dblFuncValues[1], 2.0), this->r, this->dStrike);
				}
				else if(!curPoint->isInnerPoint() && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).rightBoundary)
				{
					// Dirichlet condition when v -> inf is that U = S
					alpha[i] = tmp;
				}
				else if(!curPoint->isInnerPoint() && dblFuncValues[0] == this->myBoundingBox->getBoundary(0).rightBoundary)
				{
					// Set boundary to be the linear function at s_max
					double normalPayoff = std::max<double>(((tmp/static_cast<double>(numAssets))-strike), 0.0);
					double vRange = this->myBoundingBox->getBoundary(1).rightBoundary - this->myBoundingBox->getBoundary(1).leftBoundary;
					double sPayoffDiff = exp(dblFuncValues[0]) - normalPayoff;
					alpha[i] = normalPayoff + ((dblFuncValues[1] - this->myBoundingBox->getBoundary(1).leftBoundary) /vRange)*sPayoffDiff;
				}
				else
					alpha[i] = std::max<double>(((tmp/static_cast<double>(numAssets))-strike), 0.0);



			}
			else if (payoffType == "std_euro_put" || payoffType == "std_amer_put")
			{
				tmp = 0.0;
				for (size_t j = 0; j < numAssets; j++)
				{
					tmp += exp(dblFuncValues[2*j]);
				}
				alpha[i] = std::max<double>(strike-((tmp/static_cast<double>(numAssets))), 0.0);
			}
			else
			{
				throw new application_exception("HestonSolver::initLogTransformedGridWithPayoff : An unknown payoff-type was specified!");
			}

			delete[] dblFuncValues;
		}

		OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("HestonSolver::initLogTransformedGridWithPayoff : A grid wasn't constructed before!");
	}
}

void HestonSolver::initPATTransformedGridWithPayoff(DataVector& alpha, double strike, std::string payoffType)
{
	double tmp;

	if (this->bGridConstructed)
	{
		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
			std::stringstream coordsStream(coords);
			double* dblFuncValues = new double[dim];

			for (size_t j = 0; j < this->dim; j++)
			{
				coordsStream >> tmp;

				dblFuncValues[j] = tmp;
			}

			if (payoffType == "std_euro_call")
			{
				tmp = 0.0;
				for (size_t j = 0; j < dim; j++)
				{
					double inner_tmp = 0.0;
					for (size_t l = 0; l < dim; l++)
					{
						inner_tmp += this->eigvec_covar->get(j, l)*dblFuncValues[l];
					}
					tmp += exp(inner_tmp);
				}
				alpha[i] = std::max<double>(((tmp/static_cast<double>(dim))-strike), 0.0);
			}
			else if (payoffType == "std_euro_put" || payoffType == "std_amer_put")
			{
				tmp = 0.0;
				for (size_t j = 0; j < dim; j++)
				{
					double inner_tmp = 0.0;
					for (size_t l = 0; l < dim; l++)
					{
						inner_tmp += this->eigvec_covar->get(j, l)*dblFuncValues[l];
					}
					tmp += exp(inner_tmp);
				}
				alpha[i] = std::max<double>(strike-((tmp/static_cast<double>(dim))), 0.0);
			}
			else
			{
				throw new application_exception("HestonSolver::initPATTransformedGridWithPayoff : An unknown payoff-type was specified!");
			}

			delete[] dblFuncValues;
		}

		OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("HestonSolver::initPATTransformedGridWithPayoff : A grid wasn't constructed before!");
	}
}

double HestonSolver::evalOption(std::vector<double>& eval_point, sg::base::DataVector& alpha)
{
	std::vector<double> trans_eval = eval_point;

	// apply needed coordinate transformations
	if (this->useLogTransform)
	{
		if (this->usePAT)
		{
			for (size_t i = 0; i < eval_point.size(); i++)
			{
				double trans_point = 0.0;
				for (size_t j = 0; j < this->dim; j++)
				{
					trans_point += (this->eigvec_covar->get(j,i)*(log(eval_point[j])));
				}
				trans_point += (this->current_time*this->mu_hat->get(i));

				trans_eval[i] = trans_point;
			}
		}
		else
		{
			for (size_t i = 0; i < eval_point.size(); i=i+2)	// Here we know that the variance dimension is not log-transformed, so we shouldn't log its value
			{
				trans_eval[i] = log(trans_eval[i]);
			}
		}
	}

	sg::base::OperationEval* myEval = sg::op_factory::createOperationEval(*this->myGrid);
	double result = myEval->eval(alpha, trans_eval);
	delete myEval;

	// discounting, if PAT is used
	if (this->usePAT == true && this->payoffType != "std_amer_put")
	{
		result *= exp(((-1.0)*(this->r*this->current_time)));
	}

	return result;
}

void HestonSolver::transformPoint(sg::base::DataVector& point)
{
	sg::base::DataVector tmp_point(point);

	// apply needed coordinate transformations
	if (this->useLogTransform)
	{
		if (this->usePAT)
		{
			for (size_t i = 0; i < point.getSize(); i++)
			{
				double trans_point = 0.0;
				for (size_t j = 0; j < point.getSize(); j++)
				{
					trans_point += (this->eigvec_covar->get(j,i)*(log(point[j])));
				}
				trans_point += (this->current_time*this->mu_hat->get(i));

				tmp_point[i] = trans_point;
			}
		}
		else
		{
			for (size_t i = 0; i < point.getSize(); i++)
			{
				tmp_point[i] = log(point[i]);
			}
		}
	}

	point = tmp_point;
}

void HestonSolver::printSparseGridPAT(sg::base::DataVector& alpha, std::string tfilename, bool bSurplus) const
{
	DataVector temp(alpha);
	double tmp = 0.0;
	size_t dim = myGrid->getStorage()->dim();
	std::ofstream fileout;

	// Do Dehierarchisation, is specified
	if (bSurplus == false)
	{
		OperationHierarchisation* myHier = sg::op_factory::createOperationHierarchisation(*myGrid);
		myHier->doDehierarchisation(temp);
		delete myHier;
	}

	// Open filehandle
	fileout.open(tfilename.c_str());
	for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
	{
		std::string coords =  myGrid->getStorage()->get(i)->getCoordsStringBB(*myGrid->getBoundingBox());
		std::stringstream coordsStream(coords);

		double* dblFuncValues = new double[dim];

		for (size_t j = 0; j < dim; j++)
		{
			coordsStream >> tmp;
			dblFuncValues[j] = tmp;
		}

		for (size_t l = 0; l < dim; l++)
		{
			double trans_point = 0.0;
			for (size_t j = 0; j < dim; j++)
			{
				trans_point += this->eigvec_covar->get(l,j)*(dblFuncValues[j] - (this->current_time*this->mu_hat->get(j)));
			}
			fileout << exp(trans_point) << " ";
		}

		fileout << temp[i] << std::endl;

		delete[] dblFuncValues;
	}
	fileout.close();
}

void HestonSolver::resetSolveTime()
{
	this->current_time = 0.0;
}

size_t HestonSolver::getNeededIterationsToSolve()
{
	return this->nNeededIterations;
}

double HestonSolver::getNeededTimeToSolve()
{
	return this->dNeededTime;
}

size_t HestonSolver::getStartInnerGridSize()
{
	return this->staInnerGridSize;
}

size_t HestonSolver::getFinalInnerGridSize()
{
	return this->finInnerGridSize;
}

size_t HestonSolver::getAverageInnerGridSize()
{
	return this->avgInnerGridSize;
}

void HestonSolver::storeInnerMatrix(DataVector& alpha, std::string tFilename, double timestepsize)
{
	//    if (this->bGridConstructed)
	//      {
	//        OperationParabolicPDESolverSystemDirichlet* myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//        SGppStopwatch* myStopwatch = new SGppStopwatch();
	//
	//        std::string mtx = "";
	//
	//        myStopwatch->start();
	//        std::cout << "Generating matrix in MatrixMarket format..." << std::endl;
	//        myBSSystem->getInnerMatrix(mtx);
	//
	//        std::ofstream outfile(tFilename.c_str());
	//        outfile << mtx;
	//        outfile.close();
	//        std::cout << "Generating matrix in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;
	//
	//        delete myStopwatch;
	//        delete myBSSystem;
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::storeInnerMatrix : A grid wasn't constructed before!");
	//      }
}

void HestonSolver::storeInnerMatrixDiagonal(DataVector& alpha, std::string tFilename, double timestepsize)
{
	//    if (this->bGridConstructed)
	//      {
	//        OperationParabolicPDESolverSystemDirichlet* myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//        SGppStopwatch* myStopwatch = new SGppStopwatch();
	//
	//        std::string mtx = "";
	//
	//        myStopwatch->start();
	//        std::cout << "Generating systemmatrix's diagonal in MatrixMarket format..." << std::endl;
	//        myBSSystem->getInnerMatrixDiagonal(mtx);
	//
	//        std::ofstream outfile(tFilename.c_str());
	//        outfile << mtx;
	//        outfile.close();
	//        std::cout << "Generating systemmatrix's diagonal in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;
	//
	//        delete myStopwatch;
	//        delete myBSSystem;
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::storeInnerMatrix : A grid wasn't constructed before!");
	//      }
}

void HestonSolver::storeInnerMatrixDiagonalRowSum(DataVector& alpha, std::string tFilename, double timestepsize)
{
	//    if (this->bGridConstructed)
	//      {
	//        OperationParabolicPDESolverSystemDirichlet* myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//        SGppStopwatch* myStopwatch = new SGppStopwatch();
	//
	//        std::string mtx = "";
	//
	//        myStopwatch->start();
	//        std::cout << "Generating systemmatrix rowsum as diagonal matrix in MatrixMarket format..." << std::endl;
	//        myBSSystem->getInnerMatrixDiagonalRowSum(mtx);
	//
	//        std::ofstream outfile(tFilename.c_str());
	//        outfile << mtx;
	//        outfile.close();
	//        std::cout << "Generating systemmatrix rowsum as diagonal matrix in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;
	//
	//        delete myStopwatch;
	//        delete myBSSystem;
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::storeInnerMatrix : A grid wasn't constructed before!");
	//      }
}

void HestonSolver::storeInnerRHS(DataVector& alpha, std::string tFilename, double timestepsize)
{
	//    if (this->bGridConstructed)
	//      {
	//        OperationParabolicPDESolverSystemDirichlet* myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//        SGppStopwatch* myStopwatch = new SGppStopwatch();
	//
	//        myStopwatch->start();
	//        std::cout << "Exporting inner right-hand-side..." << std::endl;
	//        DataVector* rhs_inner = myBSSystem->generateRHS();
	//
	//        size_t nCoefs = rhs_inner->getSize();
	//        std::ofstream outfile(tFilename.c_str());
	//        for (size_t i = 0; i < nCoefs; i++)
	//          {
	//            outfile << std::scientific << rhs_inner->get(i) << std::endl;
	//          }
	//        outfile.close();
	//        std::cout << "Exporting inner right-hand-side... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;
	//
	//        delete myStopwatch;
	//        delete myBSSystem;
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::storeInnerMatrix : A grid wasn't constructed before!");
	//      }
}

void HestonSolver::storeInnerSolution(DataVector& alpha, size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, std::string tFilename)
{
	//    if (this->bGridConstructed)
	//      {
	//        Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, false, 0, myScreen);
	//        BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
	//        OperationParabolicPDESolverSystemDirichlet* myBSSystem = new BlackScholesParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
	//        SGppStopwatch* myStopwatch = new SGppStopwatch();
	//
	//        myStopwatch->start();
	//        std::cout << "Exporting inner solution..." << std::endl;
	//        myEuler->solve(*myCG, *myBSSystem, false);
	//
	//        DataVector* alpha_solve = myBSSystem->getGridCoefficientsForCG();
	//        size_t nCoefs = alpha_solve->getSize();
	//        std::ofstream outfile(tFilename.c_str());
	//        for (size_t i = 0; i < nCoefs; i++)
	//          {
	//            outfile << std::scientific << alpha_solve->get(i) << std::endl;
	//          }
	//        outfile.close();
	//
	//        std::cout << "Exporting inner solution... DONE!" << std::endl;
	//
	//        delete myStopwatch;
	//        delete myBSSystem;
	//        delete myCG;
	//        delete myEuler;
	//      }
	//    else
	//      {
	//        throw new application_exception("HestonSolver::solveImplicitEuler : A grid wasn't constructed before!");
	//      }

}

//double HestonSolver::EvaluateHestonFExact(sg::base::DataVector* characteristicFunction, double xi, double theta, double kappa, double rho, double r, double T, double K, double S, double v, int type)
//{
//	double a,b,u,x;
//	size_t vecSize = characteristicFunction->getSize();
//	//	sg::base::DataVector d(vecSize);
//	//	sg::base::DataVector g(vecSize);
//	//	sg::base::DataVector C(vecSize);
//	//	sg::base::DataVector D(vecSize);
//	//	sg::base::DataVector f(vecSize);
//	if(type == 1)
//	{
//		u = 0.5;
//		b = kappa - rho*xi;
//	}
//	else
//	{
//		u = -0.5;
//		b = kappa;
//	}
//	a = kappa*theta;
//	x = log(S);
//
//
//	double * characteristicFuncData = characteristicFunction->getPointer();
//
//	//	Build d
//	vector< complex<double> > dBuilder1;
//	vector< complex<double> > dBuilder2;
//	vector< complex<double> > d;
//	for(size_t i=0;i<vecSize;i++)
//	{
//		dBuilder1.push_back(std::complex<double>(0.0, 1.0));
//		dBuilder1[i] = pow(dBuilder1[i]*characteristicFuncData[i]*rho*xi - b,2.0);
//
//		dBuilder2.push_back(std::complex<double>(0.0, 1.0));
//		dBuilder2[i] = dBuilder2[i]*characteristicFuncData[i];
//		dBuilder2[i] = dBuilder2[i]*2.0*u - pow(characteristicFuncData[i],2.0);
//		dBuilder2[i] = dBuilder2[i]*pow(xi,2.0);
//
//		d.push_back(std::complex<double>(0.0, 0.0));
//		d[i] = sqrt(dBuilder1[i] - dBuilder2[i]);
//	}
//
//	// Build g
//	vector< complex<double> > g;
//	for(size_t i=0;i<vecSize;i++)
//	{
//		g.push_back(std::complex<double>(0.0, 1.0));
//		g[i] = (b - rho*xi*characteristicFuncData[i]*g[i] + d[i]) / (b - rho*xi*characteristicFuncData[i]*g[i] - d[i]);
//	}
//
//	// Build C
//	vector< complex<double> > C;
//	for(size_t i=0;i<vecSize;i++)
//	{
//		C.push_back(std::complex<double>(0.0, 1.0));
//		C[i] = (r*characteristicFuncData[i]*C[i]*T) + (a/pow(xi,2.0))*((b - rho*xi*characteristicFuncData[i]*C[i] + d[i])*T - 2.0*log( (1.0 - g[i]*exp(d[i]*T))/(1.0 - g[i])  ));
//	}
//
//	// Build D
//	vector< complex<double> > D;
//	for(size_t i=0;i<vecSize;i++)
//	{
//		D.push_back(std::complex<double>(0.0, 1.0));
//		D[i] = ((b - rho*xi*characteristicFuncData[i]*D[i] + d[i])/(pow(xi,2.0))) * ((1.0 - exp(d[i]*T))/(1.0 - g[i]*exp(d[i]*T)));
//	}
//
//	// Build f
//	vector< complex<double> > f;
//	for(size_t i=0;i<vecSize;i++)
//	{
//		f.push_back(std::complex<double>(0.0, 1.0));
//		f[i] = exp(C[i] + D[i]*v + f[i]*characteristicFuncData[i]*x);
//	}
//
//	// Build realArgument
//	vector< complex<double> > realArgument;
//	for(size_t i=0;i<vecSize;i++)
//	{
//		realArgument.push_back(std::complex<double>(0.0, 1.0));
//		realArgument[i] = (exp(0.0-realArgument[i]*characteristicFuncData[i]*log(K))*f[i]) / (realArgument[i]*characteristicFuncData[i]);
//	}
//
//	return real(realArgument[0]);
//}

double EvaluateHestonClosedFormIntegralFunction(double phi, double xi, double theta, double kappa, double rho, double r, double T, double K, double S, double v, int type)
{
	double a,b,u,x;
	if(type == 1)
	{
		u = 0.5;
		b = kappa - rho*xi;
	}
	else
	{
		u = -0.5;
		b = kappa;
	}
	a = kappa*theta;
	x = log(S);

	//	Build d
	complex<double> dBuilder1;
	complex<double> dBuilder2;
	complex<double> d;
	dBuilder1 = std::complex<double>(0.0, 1.0);
	dBuilder1 = pow(dBuilder1*phi*rho*xi - b,2.0);
	dBuilder2 = std::complex<double>(0.0, 1.0);
	dBuilder2 = dBuilder2*phi;
	dBuilder2 = dBuilder2*2.0*u - pow(phi,2.0);
	dBuilder2 = dBuilder2*pow(xi,2.0);
	d = sqrt(dBuilder1 - dBuilder2);

	// Build g
	complex<double> g;
	g = std::complex<double>(0.0, 1.0);
	g = (b - rho*xi*phi*g + d) / (b - rho*xi*phi*g - d);


	// Build C
	complex<double> C;
	C = std::complex<double>(0.0, 1.0);
	C = (r*phi*C*T) + (a/pow(xi,2.0))*((b - rho*xi*phi*C + d)*T - 2.0*log( (1.0 - g*exp(d*T))/(1.0 - g)  ));

	// Build D
	complex<double> D;
	D = std::complex<double>(0.0, 1.0);
	D = ((b - rho*xi*phi*D + d)/(pow(xi,2.0))) * ((1.0 - exp(d*T))/(1.0 - g*exp(d*T)));

	// Build f
	complex<double> f;
	f = complex<double>(0.0, 1.0);
	f = exp(C + D*v + f*phi*x);


	// Build realArgument
	complex<double> realArgument;
	realArgument = std::complex<double>(0.0, 1.0);
	realArgument = (exp(0.0-realArgument*phi*log(K))*f) / (realArgument*phi);

	return real(realArgument);
}

void HestonSolver::EvaluateHestonExactSurface(DataVector& alpha, double maturity)
{
	if (!this->bGridConstructed)
		throw new application_exception("HestonSolver::EvaluateHestonPriceExact : The grid wasn't initialized before!");

	if (this->numAssets != 1 || this->payoffType != "std_euro_call")
		throw new application_exception("HestonSolver::EvaluateHestonPriceExact : Can only solve in closed form for a European call option with one asset!");

	double tmp;

	for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
	{
		std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
		std::stringstream coordsStream(coords);
		double* dblFuncValues = new double[dim];

		for (size_t j = 0; j < this->dim; j++)
		{
			coordsStream >> tmp;

			if(this->useLogTransform && j == 0)
				dblFuncValues[j] = exp(tmp);
			else
				dblFuncValues[j] = tmp;
		}

		alpha[i] = EvaluateHestonPriceExact(dblFuncValues[0], dblFuncValues[1], this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0,1), this->r, maturity, this->dStrike) ;
		delete dblFuncValues;
	}

	OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
	myHierarchisation->doHierarchisation(alpha);
	delete myHierarchisation;
}

void HestonSolver::CompareHestonBs1d(double maturity, double v)
{
	if (this->numAssets != 1 || this->payoffType != "std_euro_call")
		throw new application_exception("HestonSolver::EvaluateHestonPriceExact : Can only solve in closed form for a European call option with one asset!");

	size_t dim1d = 1;
	size_t levels1d = this->levels;

	// Build a new 1d grid for stock price
	DimensionBoundary* boundaries1d = new sg::base::DimensionBoundary[dim1d];

	boundaries1d[0].leftBoundary = this->myBoundingBox->getBoundary(0).leftBoundary;
	boundaries1d[0].rightBoundary = this->myBoundingBox->getBoundary(0).rightBoundary;
	boundaries1d[0].bDirichletLeft = this->myBoundingBox->getBoundary(0).bDirichletLeft;
	boundaries1d[0].bDirichletRight = this->myBoundingBox->getBoundary(0).bDirichletRight;

	BoundingBox* boundingBox1d = new BoundingBox(dim1d, boundaries1d);

	Grid* grid1d = new LinearTrapezoidBoundaryGrid(*boundingBox1d);

	GridGenerator* myGenerator = grid1d->createGridGenerator();
	myGenerator->regular(levels1d);

	sg::base::DataVector* alphaHeston = new sg::base::DataVector(grid1d->getSize());
	sg::base::DataVector* alphaBS = new sg::base::DataVector(grid1d->getSize());

	EvaluateHestonExact1d(*alphaHeston, grid1d, boundingBox1d, maturity, v);
	EvaluateBsExact1d(*alphaBS, grid1d, boundingBox1d, maturity, 0.259);

	// Switch over the grids so we can print the 1d ones
	Grid* gridTemp = this->myGrid;
	this->myGrid = grid1d;

	printGrid(*alphaHeston, 50, "hestonExact1d.gnuplot");
	printGrid(*alphaBS, 50, "bsExact1d.gnuplot");

	alphaHeston->sub(*alphaBS);

	printGrid(*alphaHeston, 50, "hestonMinusBsExact1d.gnuplot");

	this->myGrid = gridTemp;

	delete boundaries1d;
	delete grid1d;
	delete myGenerator;
	delete alphaHeston;
	delete alphaBS;
}

/**
 * Evaluates the exact Heston solution for only varying stock price.
 */
void HestonSolver::EvaluateHestonExact1d(DataVector& alpha, Grid* grid1d, BoundingBox* boundingBox1d, double maturity, double v)
{
	double tmp;

	for (size_t i = 0; i < grid1d->getStorage()->size(); i++)
	{
		std::string coords = grid1d->getStorage()->get(i)->getCoordsStringBB(*boundingBox1d);
		std::stringstream coordsStream(coords);
		double* dblFuncValues = new double[1];
		coordsStream >> tmp;
		dblFuncValues[0] = tmp;

		alpha[i] = EvaluateHestonPriceExact(dblFuncValues[0], v, this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0,1), this->r, maturity, this->dStrike) ;
		delete dblFuncValues;
	}

	OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*grid1d);
	myHierarchisation->doHierarchisation(alpha);
	delete myHierarchisation;
}

/**
 * Evaluates the exact BS solution for only varying stock price.
 */
void HestonSolver::EvaluateBsExact1d(DataVector& alpha, Grid* grid1d, BoundingBox* boundingBox1d, double maturity, double sigma)
{
	double tmp;

	sg::finance::BlackScholesSolver* myBSSolver = new sg::finance::BlackScholesSolver(false);

	for (size_t i = 0; i < grid1d->getStorage()->size(); i++)
	{
		std::string coords = grid1d->getStorage()->get(i)->getCoordsStringBB(*boundingBox1d);
		std::stringstream coordsStream(coords);
		double* dblFuncValues = new double[1];
		coordsStream >> tmp;
		dblFuncValues[0] = tmp;

		alpha[i] = myBSSolver->getAnalyticSolution1D(dblFuncValues[0], true, maturity, sigma, this->r, this->dStrike);
		delete dblFuncValues;
	}

	OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*grid1d);
	myHierarchisation->doHierarchisation(alpha);
	delete myHierarchisation;
	delete myBSSolver;
}

double HestonSolver::EvaluateHestonPriceExact(double S, double v, double xi, double theta, double kappa, double rho, double r, double T, double K)
{
	double int1 = 0.5 + (1.0/M_PI)*GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 1);
	double int2 = 0.5 + (1.0/M_PI)*GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 2);
	return S*int1 - K*exp((-1.0)*r*T)*int2;
}

void HestonSolver::CompareHestonBsExact(sg::base::DataVector& alpha, double maturity)
{
	DataVector* alphaHeston = new DataVector(getNumberGridPoints());
	DataVector* alphaBS = new DataVector(getNumberGridPoints());

	EvaluateHestonExactSurface(*alphaHeston, maturity);
	sg::finance::BlackScholesSolver* myBSSolver = new sg::finance::BlackScholesSolver(false);

	double S,v;
	for(size_t i=0;i<this->myGridStorage->size();i++)
	{
		std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
		std::stringstream coordsStream(coords);
		coordsStream >> S;
		coordsStream >> v;

		(*alphaBS)[i] = myBSSolver->getAnalyticSolution1D(S, true, maturity, sqrt(v), this->r, this->dStrike);
	}

	OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
	myHierarchisation->doHierarchisation(*alphaBS);
	delete myHierarchisation;

	// Find the difference (heston - BS)
	for(size_t i=0;i<this->myGridStorage->size();i++)
	{
		alpha[i] = (*alphaHeston)[i] - (*alphaBS)[i];
	}

}


/** \brief Perform a single step of the Gauss-Lobatto integration

    \param f Function to integrate
    \param a Lower integration limit
    \param b Upper integration limit

    \param fa Value of function at the lower limit (used to save an
    evaluation when refinement is used)

    \param fa Value of function at the upper limit (used to save an
    evaluation when refinement is used)

    \param neval Number of evaluations made so far

    \param maxeval Maximum number of evalutions which should not be
    exceeded

    \param acc Required accuracy expressed in units of
    std::numeric_limits<double>::epsilon(). This allows less-than
    comparison by using addition and equality.
 */
double HestonSolver::GaussLobattoIntStep(
		double a, double b,
		double fa, double fb,
		size_t &neval,
		size_t maxeval,
		double acc
		, double xi, double theta, double kappa, double rho, double r, double T, double K, double S, double v, int type)
{

	// Constants used in the algorithm
	const double alpha = std::sqrt(2.0/3.0);
	const double beta  = 1.0/std::sqrt(5.0);

	if (neval >= maxeval)
	{
		throw new application_exception("HestonSolver::Gauss-Lobatto : Maximum number of evaluations reached in GaussLobatto.");
	}

	// Here the abcissa points and function values for both the 4-point
	// and the 7-point rule are calculated (the points at the end of
	// interval come from the function call, i.e., fa and fb. Also note
	// the 7-point rule re-uses all the points of the 4-point rule.)
	const double h=(b-a)/2;
	const double m=(a+b)/2;

	const double mll=m-alpha*h;
	const double ml =m-beta*h;
	const double mr =m+beta*h;
	const double mrr=m+alpha*h;

	const double fmll= EvaluateHestonClosedFormIntegralFunction(mll, xi, theta, kappa, rho, r, T, K, S, v, type);
	const double fml = EvaluateHestonClosedFormIntegralFunction(ml, xi, theta, kappa, rho, r, T, K, S, v, type);
	const double fm  = EvaluateHestonClosedFormIntegralFunction(m, xi, theta, kappa, rho, r, T, K, S, v, type);
	const double fmr = EvaluateHestonClosedFormIntegralFunction(mr, xi, theta, kappa, rho, r, T, K, S, v, type);
	const double fmrr= EvaluateHestonClosedFormIntegralFunction(mrr, xi, theta, kappa, rho, r, T, K, S, v, type);
	neval+=5;

	// Both the 4-point and 7-point rule integrals are evaluted
	const double integral2=(h/6)*(fa+fb+5*(fml+fmr));
	const double integral1=(h/1470)*(77*(fa+fb)
			+432*(fmll+fmrr)+625*(fml+fmr)+672*fm);

	// The difference betwen the 4-point and 7-point integrals is the
	// estimate of the accuracy
	const double estacc=(integral1-integral2);

	// The volatile keyword should prevent the floating point
	// destination value from being stored in extended precision
	// registers which actually have a very different
	// std::numeric_limits<double>::epsilon().
	volatile double dist = acc + estacc;

	if(dist==acc || mll<=a || b<=mrr)
	{
		if (not (m>a && b>m))
		{
			throw new application_exception("HestonSolver::Gauss-Lobatto : Integration reached an interval with no more machine numbers!");
		}
		return integral1;
	}
	else {
		return  GaussLobattoIntStep(a, mll, fa, fmll, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
				+ GaussLobattoIntStep(mll, ml, fmll, fml, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
				+ GaussLobattoIntStep(ml, m, fml, fm, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
				+ GaussLobattoIntStep(m, mr, fm, fmr, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
				+ GaussLobattoIntStep(mr, mrr, fmr, fmrr, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
				+ GaussLobattoIntStep(mrr, b, fmrr, fb, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type);

	}
}

/** \brief Compute the Gauss-Lobatto integral

    \param f The function to be integrated

    \param a The lower integration limit

    \param b The upper integration limit

    \param abstol Absolute tolerance -- integration stops when the
    error estimate is smaller than this

    \param maxeval Maxium of evaluations to make. If this number of
    evalution is made without reaching the requied accuracy, an
    exception of type std::runtime_error is thrown.
 */
double HestonSolver::GaussLobattoInt(double a, double b,
		double abstol,
		size_t maxeval
		, double xi, double theta, double kappa, double rho, double r, double T, double K, double S, double v, int type)
{
	const double tol_epsunit=abstol/std::numeric_limits<double>::epsilon();
	size_t neval=0;
	return GaussLobattoIntStep(a, b,
			EvaluateHestonClosedFormIntegralFunction(a, xi, theta, kappa, rho, r, T, K, S, v, type), EvaluateHestonClosedFormIntegralFunction(b, xi, theta, kappa, rho, r, T, K, S, v, type),
			neval,
			maxeval,
			tol_epsunit, xi, theta, kappa, rho, r, T, K, S, v, type);
}

void HestonSolver::CompareHestonSolutionToExact(sg::base::DataVector* solution, sg::base::DataVector* exact, std::string filename, size_t PointsPerDimension)
{
	DimensionBoundary dimOne;
	DimensionBoundary dimTwo;
	std::ofstream fileout;

	fileout.open(filename.c_str());
	OperationEval* myEval = sg::op_factory::createOperationEval(*myGrid);


	if (myGrid->getStorage()->dim() == 2)
	{
		dimOne = myGrid->getBoundingBox()->getBoundary(0);
		dimTwo = myGrid->getBoundingBox()->getBoundary(1);

		double offset_x = dimOne.leftBoundary;
		double offset_y = dimTwo.leftBoundary;
		double inc_x = ((dimOne.rightBoundary - dimOne.leftBoundary)/(PointsPerDimension-1));
		double inc_y = ((dimTwo.rightBoundary - dimTwo.leftBoundary)/(PointsPerDimension-1));

		size_t points = (size_t)PointsPerDimension;

		for (size_t i = 0; i < points; i++)
		{
			for (size_t j = 0; j < points; j++)
			{
				std::vector<double> point;
				point.push_back(offset_x+(((double)(i))*inc_x));
				point.push_back(offset_y+(((double)(j))*inc_y));
				fileout << (offset_x+((double)(i))*inc_x) << " " << (offset_y+((double)(j))*inc_y) << " " << (myEval->eval(*solution,point) - myEval->eval(*exact,point)) << std::endl;
			}
			fileout << std::endl;
		}
	}

	delete myEval;
	// close filehandle
	fileout.close();
}



}
}

