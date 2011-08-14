/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/MPI/SGppMPITools.hpp"
#include "solver/sle/BiCGStabMPI.hpp"
#include "solver/sle/ConjugateGradientsMPI.hpp"
#include "algorithm/pde/BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI.hpp"
#include "algorithm/pde/BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI.hpp"
#include "application/pde/BlackScholesSolverMPI.hpp"

#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/ode/AdamsBashforth.hpp"
#include "solver/ode/VarTimestep.hpp"
#include "solver/ode/StepsizeControlH.hpp"
#include "solver/ode/StepsizeControlBDF.hpp"
#include "solver/ode/StepsizeControlEJ.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include "basis/operations_factory.hpp"
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace sg
{
namespace parallel
{

BlackScholesSolverMPI::BlackScholesSolverMPI(bool useLogTransform, bool usePAT) : sg::pde::ParabolicPDESolver()
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

BlackScholesSolverMPI::~BlackScholesSolverMPI()
{
	if (this->bStochasticDataAlloc)
	{
		delete this->mus;
		delete this->sigmas;
		delete this->rhos;
	}
	if (this->myScreen != NULL)
	{
		delete this->myScreen;
	}
}

void BlackScholesSolverMPI::getGridNormalDistribution(sg::base::DataVector& alpha, std::vector<double>& norm_mu, std::vector<double>& norm_sigma)
{
	if (this->bGridConstructed)
	{
		double tmp;
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
		throw new application_exception("BlackScholesSolverMPI::getGridNormalDistribution : The grid wasn't initialized before!");
	}
}

void BlackScholesSolverMPI::constructGrid(sg::base::BoundingBox& BoundingBox, size_t level)
{
	this->dim = BoundingBox.getDimensions();
	this->levels = level;

	this->myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

	sg::base::GridGenerator* myGenerator = this->myGrid->createGridGenerator();
	myGenerator->regular(this->levels);
	delete myGenerator;

	this->myBoundingBox = this->myGrid->getBoundingBox();
	this->myGridStorage = this->myGrid->getStorage();

	//std::string serGrid;
	//myGrid->serialize(serGrid);
	//std::cout << serGrid << std::endl;

	this->bGridConstructed = true;
}

void BlackScholesSolverMPI::setStochasticData(sg::base::DataVector& mus, sg::base::DataVector& sigmas, sg::base::DataMatrix& rhos, double r)
{
	this->mus = new sg::base::DataVector(mus);
	this->sigmas = new sg::base::DataVector(sigmas);
	this->rhos = new sg::base::DataMatrix(rhos);
	this->r = r;

	// calculate eigenvalues, eigenvectors and mu_hat from stochastic data for PAT
	size_t mydim = this->mus->getSize();
	this->eigval_covar = new sg::base::DataVector(mydim);
	this->eigvec_covar = new sg::base::DataMatrix(mydim,mydim);
	this->mu_hat = new sg::base::DataVector(mydim);

	// 1d test case
	if (mydim == 1)
	{
		this->eigval_covar->set(0, this->sigmas->get(0)*this->sigmas->get(0));
		this->eigvec_covar->set(0, 0, 1.0);
	}
	// 2d test case
	if (mydim == 2)
	{
		// correlation -0.5, sigma_1 0.3, sigma_2 0.4
		this->eigval_covar->set(0, 0.0555377800527509792);
		this->eigval_covar->set(1, 0.194462219947249021);

		this->eigvec_covar->set(0, 0, -0.867142152569025494);
		this->eigvec_covar->set(0, 1, 0.498060726456078796);
		this->eigvec_covar->set(1, 0, -0.498060726456078796);
		this->eigvec_covar->set(1, 1, -0.867142152569025494);

		// correlation 0.1, sigma 0.4
//		this->eigval_covar->set(0, 0.176);
//		this->eigval_covar->set(1, 0.144);
//
//		this->eigvec_covar->set(0, 0, 0.707106781186547351);
//		this->eigvec_covar->set(0, 1, -0.707106781186547573);
//		this->eigvec_covar->set(1, 0, 0.707106781186547573);
//		this->eigvec_covar->set(1, 1, 0.707106781186547351);

		// correlation 0.25, sigma 0.4
//		this->eigval_covar->set(0, 0.20);
//		this->eigval_covar->set(1, 0.12);
//
//		this->eigvec_covar->set(0, 0, 0.707106781186547573);
//		this->eigvec_covar->set(0, 1, -0.707106781186547462);
//		this->eigvec_covar->set(1, 0, 0.707106781186547462);
//		this->eigvec_covar->set(1, 1, 0.707106781186547573);

		// correlation 0.5, sigma 0.4
//		this->eigval_covar->set(0, 0.24);
//		this->eigval_covar->set(1, 0.08);
//
//		this->eigvec_covar->set(0, 0, 0.707106781186547462);
//		this->eigvec_covar->set(0, 1, -0.707106781186547462);
//		this->eigvec_covar->set(1, 0, 0.707106781186547462);
//		this->eigvec_covar->set(1, 1, 0.707106781186547462);

		// correlation -0.5, sigma 0.4
//		this->eigval_covar->set(0, 0.24);
//		this->eigval_covar->set(1, 0.08);
//
//		this->eigvec_covar->set(0, 0, 0.707106781186547462);
//		this->eigvec_covar->set(0, 1, 0.707106781186547462);
//		this->eigvec_covar->set(1, 0, -0.707106781186547462);
//		this->eigvec_covar->set(1, 1, 0.707106781186547462);

		// correlation 0.0, sigma 0.4
//		this->eigval_covar->set(0, 0.16);
//		this->eigval_covar->set(1, 0.16);
//
//		this->eigvec_covar->set(0, 0, 1);
//		this->eigvec_covar->set(0, 1, 0);
//		this->eigvec_covar->set(1, 0, 0);
//		this->eigvec_covar->set(1, 1, 1);
	}
	// 3d test case
	if (mydim == 3)
	{
		this->eigval_covar->set(0, 0.0161152062670340546);
		this->eigval_covar->set(1, 0.109759129882028184);
		this->eigval_covar->set(2, 0.164125663850937770);

		this->eigvec_covar->set(0, 0, -0.869836297894464927);
		this->eigvec_covar->set(0, 1, -0.472520156399458380);
		this->eigvec_covar->set(0, 2, -0.141808027493095845);
		this->eigvec_covar->set(1, 0, -0.493287590771869233);
		this->eigvec_covar->set(1, 1, 0.837245861704536964);
		this->eigvec_covar->set(1, 2, 0.235980337844304722);
		this->eigvec_covar->set(2, 0, -0.722271802968979613e-2);
		this->eigvec_covar->set(2, 1, -0.275216403680555832);
		this->eigvec_covar->set(2, 2, 0.961355170313971441);
	}
	// 4d test case
	if (mydim == 4)
	{
		this->eigval_covar->set(0, 0.203896808612126890);
		this->eigval_covar->set(1, 0.143228838600683389);
		this->eigval_covar->set(2, 0.0369289052706551282);
		this->eigval_covar->set(3, 0.884454475165345477e-1);

		this->eigvec_covar->set(0, 0, -0.758784156527507303);
		this->eigvec_covar->set(0, 1, -0.441609644035335480);
		this->eigvec_covar->set(0, 2, -0.329723883556484409);
		this->eigvec_covar->set(0, 3, 0.347145051398191185);
		this->eigvec_covar->set(1, 0, 0.0381555338704252095);
		this->eigvec_covar->set(1, 1, -0.502978780370165746e-1);
		this->eigvec_covar->set(1, 2, 0.713853106787664671);
		this->eigvec_covar->set(1, 3, 0.697443919343796126);
		this->eigvec_covar->set(2, 0, 0.327315044904916141);
		this->eigvec_covar->set(2, 1, 0.376967846478724333);
		this->eigvec_covar->set(2, 2, -0.600862700798999505);
		this->eigvec_covar->set(2, 3, 0.624278879098610684);
		this->eigvec_covar->set(3, 0, -0.561832377508447611);
		this->eigvec_covar->set(3, 1, 0.812616938342507478);
		this->eigvec_covar->set(3, 2, 0.143735581296216858);
		this->eigvec_covar->set(3, 3, -0.577769311359955864e-1);
	}
	// 5d test case
	if (mydim == 5)
	{
		this->eigval_covar->set(0, 0.248090694157781677);
		this->eigval_covar->set(1, 0.181003240820949346);
		this->eigval_covar->set(2, 0.0132179416451147155);
		this->eigval_covar->set(3, 0.0939549786605669707);
		this->eigval_covar->set(4, 0.0587331447155870698);

		this->eigvec_covar->set(0, 0, 0.263790550378285305);
		this->eigvec_covar->set(0, 1, 0.859923642135395405);
		this->eigvec_covar->set(0, 2, 0.334323851170300113);
		this->eigvec_covar->set(0, 3, -0.280069559202008711);
		this->eigvec_covar->set(0, 4, 0.0271012873267933163);
		this->eigvec_covar->set(1, 0, -0.460592496016481306e-1);
		this->eigvec_covar->set(1, 1, 0.127520536615915014e-1);
		this->eigvec_covar->set(1, 2, -0.377657780213884409);
		this->eigvec_covar->set(1, 3, -0.528472957191017390);
		this->eigvec_covar->set(1, 4, -0.758819389061222926);
		this->eigvec_covar->set(2, 0, -0.164069200805689486);
		this->eigvec_covar->set(2, 1, -0.378822519988808670);
		this->eigvec_covar->set(2, 2, 0.470800766942319815);
		this->eigvec_covar->set(2, 3, -0.728869775137000020);
		this->eigvec_covar->set(2, 4, 0.276893994941337263);
		this->eigvec_covar->set(3, 0, 0.749451659161123107);
		this->eigvec_covar->set(3, 1, -0.143542996441549164);
		this->eigvec_covar->set(3, 2, -0.467549699477188219);
		this->eigvec_covar->set(3, 3, -0.257715994147244276);
		this->eigvec_covar->set(3, 4, 0.364276493384776856);
		this->eigvec_covar->set(4, 0, -0.582834967194725273);
		this->eigvec_covar->set(4, 1, 0.310254123817747751);
		this->eigvec_covar->set(4, 2, -0.552581288090625233);
		this->eigvec_covar->set(4, 3, -0.211207700566551748);
		this->eigvec_covar->set(4, 4, 0.462699694124277694);
	}

	for (size_t i = 0; i < mydim; i++)
	{
		double tmp = 0.0;
		for (size_t j = 0; j < mydim; j++)
		{
			tmp += ((this->mus->get(j) - (0.5*this->sigmas->get(j)*this->sigmas->get(j))) * this->eigvec_covar->get(j, i));
		}
		this->mu_hat->set(i, tmp);
	}

	bStochasticDataAlloc = true;
}

void BlackScholesSolverMPI::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed && this->bStochasticDataAlloc)
	{
		solver::Euler* myEuler = new solver::Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
		solver::SLESolver* myCG;
		pde::OperationParabolicPDESolverSystem* myBSSystem = NULL;

		if (this->usePAT == false)
		{
			myCG = new parallel::BiCGStabMPI(maxCGIterations, epsilonCG);
			myBSSystem = new parallel::BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ExEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
		}
		else
		{
			myCG = new parallel::ConjugateGradientsMPI(maxCGIterations, epsilonCG);
			myBSSystem = new parallel::BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI(*this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, timestepsize, "ExEul", this->dStrike, this->payoffType, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
		}

		base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
		this->staInnerGridSize = getNumberInnerGridPoints();

		myStopwatch->start();
		myEuler->solve(*myCG, *myBSSystem, true, verbose);
		this->dNeededTime = myStopwatch->stop();

		if (myGlobalMPIComm->getMyRank() == 0)
		{
			std::cout << "Using Explicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
			std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
			std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;
			std::cout << "Average Grid size: " << static_cast<double>(myBSSystem->getSumGridPointsComplete())/static_cast<double>(numTimesteps) << std::endl;
			std::cout << "Average Grid size (Inner): " << static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

			if (this->myScreen != NULL)
			{
				std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
				this->myScreen->writeEmptyLines(2);
			}
		}

		this->finInnerGridSize = getNumberInnerGridPoints();
		this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps))+0.5);
		this->nNeededIterations = myEuler->getNumberIterations();

		delete myBSSystem;
		delete myCG;
		delete myEuler;
		delete myStopwatch;

		this->current_time += (static_cast<double>(numTimesteps)*timestepsize);
	}
	else
	{
		throw new application_exception("BlackScholesSolverMPI::solveExplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
	}
}

void BlackScholesSolverMPI::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed && this->bStochasticDataAlloc)
	{
		solver::Euler* myEuler = new solver::Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
		solver::SLESolver* myCG;
		pde::OperationParabolicPDESolverSystem* myBSSystem = NULL;

		if (this->usePAT == false)
		{
			myCG = new parallel::BiCGStabMPI(maxCGIterations, epsilonCG);
			myBSSystem = new parallel::BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
		}
		else
		{
			myCG = new parallel::ConjugateGradientsMPI(maxCGIterations, epsilonCG);
			myBSSystem = new parallel::BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI(*this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
		}

		base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
		this->staInnerGridSize = getNumberInnerGridPoints();

		myStopwatch->start();
		myEuler->solve(*myCG, *myBSSystem, true, verbose);
		this->dNeededTime = myStopwatch->stop();

		if (myGlobalMPIComm->getMyRank() == 0)
		{
			std::cout << "Using Implicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
			std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
			std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;
			std::cout << "Average Grid size: " << static_cast<double>(myBSSystem->getSumGridPointsComplete())/static_cast<double>(numTimesteps) << std::endl;
			std::cout << "Average Grid size (Inner): " << static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

			if (this->myScreen != NULL)
			{
				std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
				this->myScreen->writeEmptyLines(2);
			}
		}

		this->finInnerGridSize = getNumberInnerGridPoints();
		this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps))+0.5);
		this->nNeededIterations = myEuler->getNumberIterations();

		delete myBSSystem;
		delete myCG;
		delete myEuler;
		delete myStopwatch;

		this->current_time += (static_cast<double>(numTimesteps)*timestepsize);
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solveImplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
	}
}

void BlackScholesSolverMPI::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, size_t NumImEul)
{
	if (this->bGridConstructed && this->bStochasticDataAlloc)
	{
		solver::SLESolver* myCG;
		pde::OperationParabolicPDESolverSystem* myBSSystem = NULL;

		if (this->usePAT == false)
		{
			myCG = new parallel::BiCGStabMPI(maxCGIterations, epsilonCG);
			myBSSystem = new parallel::BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
		}
		else
		{
			myCG = new parallel::ConjugateGradientsMPI(maxCGIterations, epsilonCG);
			myBSSystem = new parallel::BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI(*this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, timestepsize, "ImEul", this->dStrike, this->payoffType, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
		}

		base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
		this->staInnerGridSize = getNumberInnerGridPoints();

		size_t numCNSteps;
		size_t numIESteps;

		numCNSteps = numTimesteps;
		if (numTimesteps > NumImEul)
		{
			numCNSteps = numTimesteps - NumImEul;
		}
		numIESteps = NumImEul;

		solver::Euler* myEuler = new solver::Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
		solver::CrankNicolson* myCN = new solver::CrankNicolson(numCNSteps, timestepsize, this->myScreen);

		myStopwatch->start();
		if (numIESteps > 0)
		{
			if (myGlobalMPIComm->getMyRank() == 0)
			{
				std::cout << "Using Implicit Euler to solve " << numIESteps << " timesteps:" << std::endl;
			}
			myBSSystem->setODESolver("ImEul");
			myEuler->solve(*myCG, *myBSSystem, false, false);
		}

		if (myGlobalMPIComm->getMyRank() == 0)
		{
			std::cout << "Using Crank Nicolson to solve " << numCNSteps << " timesteps:" << std::endl << std::endl << std::endl << std::endl;
		}

		myBSSystem->setODESolver("CrNic");
		myCN->solve(*myCG, *myBSSystem, true, false);
		this->dNeededTime = myStopwatch->stop();

		if (myGlobalMPIComm->getMyRank() == 0)
		{
			std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
			std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;
			std::cout << "Average Grid size: " << static_cast<double>(myBSSystem->getSumGridPointsComplete())/static_cast<double>(numTimesteps) << std::endl;
			std::cout << "Average Grid size (Inner): " << static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

			if (this->myScreen != NULL)
			{
				std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
				this->myScreen->writeEmptyLines(2);
			}
		}

		this->finInnerGridSize = getNumberInnerGridPoints();
		this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myBSSystem->getSumGridPointsInner())/static_cast<double>(numTimesteps))+0.5);
		this->nNeededIterations = myEuler->getNumberIterations() + myCN->getNumberIterations();

		delete myBSSystem;
		delete myCG;
		delete myCN;
		delete myEuler;
		delete myStopwatch;

		this->current_time += (static_cast<double>(numTimesteps)*timestepsize);
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solveCrankNicolson : A grid wasn't constructed before or stochastic parameters weren't set!");
	}
}


void BlackScholesSolverMPI::solveAdamsBashforth(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose)
{
	throw new application_exception("BlackScholesSolverMPI::solveAdamsBashforth : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveSCAC(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose)
{
	throw new application_exception("BlackScholesSolverMPI::solveSCAC : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveSCH(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose)
{
	throw new application_exception("BlackScholesSolverMPI::solveSCH : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveSCBDF(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose)
{
	throw new application_exception("BlackScholesSolverMPI::solveSCBDF : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveSCEJ(size_t numTimesteps, double timestepsize, double epsilon, double myAlpha, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose)
{
	throw new application_exception("BlackScholesSolverMPI::solveSCEJ : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveX(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose, void *myODESolverV, std::string Solver)
{
	throw new application_exception("BlackScholesSolverMPI::solveX : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::initGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType)
{
	this->dStrike = strike;
	this->payoffType = payoffType;

	if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put")
	{
		this->tBoundaryType = "Dirichlet";
	}

	if (this->useLogTransform)
	{
		if (this->usePAT == true)
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
	}
}

double BlackScholesSolverMPI::get1DEuroCallPayoffValue(double assetValue, double strike)
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

std::vector<size_t> BlackScholesSolverMPI::getAlgorithmicDimensions()
{
	return this->myGrid->getAlgorithmicDimensions();
}

void BlackScholesSolverMPI::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims)
{
	if (this->tBoundaryType == "freeBoundaries")
	{
		this->myGrid->setAlgorithmicDimensions(newAlgoDims);
	}
	else
	{
		throw new application_exception("BlackScholesSolver::setAlgorithmicDimensions : Set algorithmic dimensions is only supported when choosing option type all!");
	}
}

void BlackScholesSolverMPI::initScreen()
{
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Black Scholes Solver, 2.0.0", "TUM (C) 2009-2010, by Alexander Heinecke");
	this->myScreen->writeStartSolve("Multidimensional Black Scholes Solver");
}

void BlackScholesSolverMPI::setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, size_t refineMaxLevel, int numCoarsenPoints, double coarsenThreshold, double refineThreshold)
{
	this->useCoarsen = true;
	this->coarsenThreshold = coarsenThreshold;
	this->refineThreshold = refineThreshold;
	this->refineMaxLevel = refineMaxLevel;
	this->adaptSolveMode = adaptSolveMode;
	this->refineMode = refineMode;
	this->numCoarsenPoints = numCoarsenPoints;
}

void BlackScholesSolverMPI::initCartesianGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType)
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
					tmp += dblFuncValues[j];
				}
				alpha[i] = std::max<double>(((tmp/static_cast<double>(dim))-strike), 0.0);
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
				throw new application_exception("BlackScholesSolver::initCartesianGridWithPayoff : An unknown payoff-type was specified!");
			}

			delete[] dblFuncValues;
		}

		base::OperationHierarchisation* myHierarchisation = sg::GridOperationFactory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::initCartesianGridWithPayoff : A grid wasn't constructed before!");
	}
}

void BlackScholesSolverMPI::initLogTransformedGridWithPayoff(DataVector& alpha, double strike, std::string payoffType)
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
					tmp += exp(dblFuncValues[j]);
				}
				alpha[i] = std::max<double>(((tmp/static_cast<double>(dim))-strike), 0.0);
			}
			else if (payoffType == "std_euro_put" || payoffType == "std_amer_put")
			{
				tmp = 0.0;
				for (size_t j = 0; j < dim; j++)
				{
					tmp += exp(dblFuncValues[j]);
				}
				alpha[i] = std::max<double>(strike-((tmp/static_cast<double>(dim))), 0.0);
			}
			else
			{
				throw new application_exception("BlackScholesSolver::initLogTransformedGridWithPayoff : An unknown payoff-type was specified!");
			}

			delete[] dblFuncValues;
		}

		base::OperationHierarchisation* myHierarchisation = sg::GridOperationFactory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::initLogTransformedGridWithPayoff : A grid wasn't constructed before!");
	}
}

void BlackScholesSolverMPI::initPATTransformedGridWithPayoff(DataVector& alpha, double strike, std::string payoffType)
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
				throw new application_exception("BlackScholesSolverMPI::initPATTransformedGridWithPayoff : An unknown payoff-type was specified!");
			}

			delete[] dblFuncValues;
		}

		OperationHierarchisation* myHierarchisation = sg::GridOperationFactory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("BlackScholesSolverMPI::initPATTransformedGridWithPayoff : A grid wasn't constructed before!");
	}
}

double BlackScholesSolverMPI::evalOption(std::vector<double>& eval_point, sg::base::DataVector& alpha)
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
			for (size_t i = 0; i < eval_point.size(); i++)
			{
				trans_eval[i] = log(trans_eval[i]);
			}
		}
	}

	sg::base::OperationEval* myEval = sg::GridOperationFactory::createOperationEval(*this->myGrid);
	double result = myEval->eval(alpha, trans_eval);
	delete myEval;

	// discounting, if PAT is used
	if (this->usePAT == true)
	{
		result *= exp(((-1.0)*(this->r*this->current_time)));
	}

	return result;
}

void BlackScholesSolverMPI::transformPoint(sg::base::DataVector& point)
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

void BlackScholesSolverMPI::printSparseGridPAT(sg::base::DataVector& alpha, std::string tfilename, bool bSurplus) const
{
	DataVector temp(alpha);
	double tmp = 0.0;
	size_t dim = myGrid->getStorage()->dim();
	std::ofstream fileout;

	// Do Dehierarchisation, is specified
	if (bSurplus == false)
	{
		OperationHierarchisation* myHier = sg::GridOperationFactory::createOperationHierarchisation(*myGrid);
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

void BlackScholesSolverMPI::resetSolveTime()
{
	this->current_time = 0.0;
}


size_t BlackScholesSolverMPI::getNeededIterationsToSolve()
{
	return this->nNeededIterations;
}

double BlackScholesSolverMPI::getNeededTimeToSolve()
{
	return this->dNeededTime;
}

size_t BlackScholesSolverMPI::getStartInnerGridSize()
{
	return this->staInnerGridSize;
}

size_t BlackScholesSolverMPI::getFinalInnerGridSize()
{
	return this->finInnerGridSize;
}

size_t BlackScholesSolverMPI::getAverageInnerGridSize()
{
	return this->avgInnerGridSize;
}

}

}
