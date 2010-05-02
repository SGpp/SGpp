/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "algorithm/pde/BlackScholesODESolverSystem.hpp"
#include "application/pde/BlackScholesSolver.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include <cstdlib>
#include <sstream>
#include <cmath>

namespace sg
{

BlackScholesSolver::BlackScholesSolver()
{
	bStochasticDataAlloc = false;
	bGridConstructed = false;
	myScreen = NULL;
}

BlackScholesSolver::~BlackScholesSolver()
{
	if (bStochasticDataAlloc)
	{
		delete mus;
		delete sigmas;
		delete rhos;
	}
	if (bGridConstructed)
	{
		delete myGrid;
	}
	if (myScreen != NULL)
	{
		delete myScreen;
	}
}

void BlackScholesSolver::constructGrid(BoundingBox& BoundingBox, size_t level)
{
	dim = BoundingBox.getDimensions();
	levels = level;

	myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

	GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(levels);
	delete myGenerator;

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	//std::string serGrid;
	//myGrid->serialize(serGrid);
	//std::cout << serGrid << std::endl;

	bGridConstructed = true;
}

void BlackScholesSolver::refineInitialGrid(DataVector& alpha, double* strike, std::string payoffType, double dStrikeDistance)
{
	size_t nRefinements = 0;

	if (bGridConstructed)
	{

		DataVector refineVector(alpha.getSize());

		if (payoffType == "avgM")
		{
			double tmp;
			double* dblFuncValues = new double[dim];
			double dDistance;

			for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
			{
				std::string coords = myGridStorage->get(i)->getCoordsStringBB(*myBoundingBox);
				std::stringstream coordsStream(coords);

				for (size_t j = 0; j < dim; j++)
				{
					coordsStream >> tmp;

					dblFuncValues[j] = tmp;
				}

				tmp = 0.0;
				for (size_t j = 0; j < dim; j++)
				{
					tmp += dblFuncValues[j];
				}
				dDistance = fabs(((tmp/static_cast<double>(dim))-1.0));
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

			myGrid->createGridGenerator()->refine(myRefineFunc);

			delete myRefineFunc;

			alpha.resize(myGridStorage->size());

			// reinit the grid with the payoff function
			initGridWithEuroCallPayoff(alpha, strike, payoffType);
		}
		else
		{
			throw new application_exception("BlackScholesSolver::refineGrid : An unsupported payoffType was specified!");
		}
	}
	else
	{
		throw new application_exception("BlackScholesSolver::refineInitialGrid : The grid wasn't initialized before!");
	}
}

void BlackScholesSolver::refineInitialGridSurplus(DataVector& alpha, double dPercentage)
{
	size_t nRefinements = 0;

	if (bGridConstructed)
	{
		SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&alpha, nRefinements, 0.0);

		myGrid->createGridGenerator()->refine(myRefineFunc);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	}
	else
	{
		throw new application_exception("BlackScholesSolver::refineIntialGridSurplus : The grid wasn't initialized before!");
	}
}

void BlackScholesSolver::constructGrid(std::string tfilename, DataVector& emptyAlpha, bool& ishierarchized)
{
	IOToolBonnSG* myImporter = new IOToolBonnSG();
	std::string serGrid;

	// test if emptyAlpha is really empty
	if (emptyAlpha.getSize() != 0)
	{
		throw new application_exception("BlackScholesSolver::constructGrid : A non-empty coefficients' vector was used in gird construction from file!");
	}

	myImporter->readFile(tfilename, serGrid, emptyAlpha, ishierarchized);

	myGrid = Grid::unserialize(serGrid);

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	// Set every boundary to dirichlet boundaries
	for (size_t i = 0; i < myGridStorage->dim(); i++)
	{
		DimensionBoundary myDimBound = myBoundingBox->getBoundary(i);

		myDimBound.bDirichletLeft = true;
		myDimBound.bDirichletRight = true;

		myBoundingBox->setBoundary(i, myDimBound);
	}

	if (ishierarchized == false)
	{
		OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
		myHierarchisation->doHierarchisation(emptyAlpha);
		delete myHierarchisation;
	}

	bGridConstructed = true;

	delete myImporter;
}

void BlackScholesSolver::storeGrid(std::string tfilename, DataVector& alpha, bool ishierarchized)
{
	IOToolBonnSG* myExporter = new IOToolBonnSG();
	DataVector copyAlpha(alpha);

	if (bGridConstructed)
	{
		if (ishierarchized == false)
		{
			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doDehierarchisation(copyAlpha);
			delete myHierarchisation;
		}

		myExporter->writeFile(tfilename, *myGrid, copyAlpha, ishierarchized);
	}
	else
	{
		throw new application_exception("BlackScholesSolver::storeGrid : A grid wasn't constructed before!");
	}

	delete myExporter;
}

void BlackScholesSolver::setStochasticData(DataVector& mus, DataVector& sigmas, DataVector& rhos, double r)
{
	this->mus = new DataVector(mus);
	this->sigmas = new DataVector(sigmas);
	this->rhos = new DataVector(rhos);
	this->r = r;

	bStochasticDataAlloc = true;
}

void BlackScholesSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (bGridConstructed && bStochasticDataAlloc)
	{
		Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
		BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
		BlackScholesODESolverSystem* myBSSystem = new BlackScholesODESolverSystem(*myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, r, timestepsize, "ExEul");
		SGppStopwatch* myStopwatch = new SGppStopwatch();
		double execTime;

		myStopwatch->start();
		myEuler->solve(*myCG, *myBSSystem, verbose);
		execTime = myStopwatch->stop();

		if (myScreen != NULL)
		{
			std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
			myScreen->writeEmptyLines(2);
		}

		delete myBSSystem;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solveExplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
	}
}

void BlackScholesSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (bGridConstructed && bStochasticDataAlloc)
	{
		Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
		BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
		BlackScholesODESolverSystem* myBSSolver = new BlackScholesODESolverSystem(*myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, r, timestepsize, "ImEul");
		SGppStopwatch* myStopwatch = new SGppStopwatch();
		double execTime;

		myStopwatch->start();
		myEuler->solve(*myCG, *myBSSolver, verbose);
		execTime = myStopwatch->stop();

		if (myScreen != NULL)
		{
			std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
			myScreen->writeEmptyLines(2);
		}

		delete myBSSolver;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solveImplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
	}
}

void BlackScholesSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha)
{
	if (bGridConstructed && bStochasticDataAlloc)
	{
		CrankNicolson* myCN = new CrankNicolson(numTimesteps, timestepsize, myScreen);
		BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
		BlackScholesODESolverSystem* myBSSolver = new BlackScholesODESolverSystem(*myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, r, timestepsize, "CrNic");
		SGppStopwatch* myStopwatch = new SGppStopwatch();
		double execTime;

		myStopwatch->start();
		myCN->solve(*myCG, *myBSSolver, false);
		execTime = myStopwatch->stop();

		if (myScreen != NULL)
		{
			std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
			myScreen->writeEmptyLines(2);
		}

		delete myBSSolver;
		delete myCG;
		delete myCN;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solveCrankNicolson : A grid wasn't constructed before or stochastic parameters weren't set!");
	}
}

void BlackScholesSolver::printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename)
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void BlackScholesSolver::printSparseGrid(DataVector& alpha, std::string tfilename, bool bSurplus)
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

void BlackScholesSolver::initGridWithEuroCallPayoff(DataVector& alpha, double* strike, std::string payoffType)
{
	double tmp;

	if (bGridConstructed)
	{
		for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
		{
			std::string coords = myGridStorage->get(i)->getCoordsStringBB(*myBoundingBox);
			std::stringstream coordsStream(coords);
			double* dblFuncValues = new double[dim];

			for (size_t j = 0; j < dim; j++)
			{
				coordsStream >> tmp;
				if ((payoffType == "max") || (payoffType == "avg"))
				{
					dblFuncValues[j] = get1DEuroCallPayoffValue(tmp, strike[j]);
				}
				if (payoffType == "avgM")
				{
					dblFuncValues[j] = tmp;
				}
			}

			if (payoffType == "max")
			{
				tmp = 0.0;
				for (size_t j = 0; j < dim; j++)
				{
					if (dblFuncValues[j] > tmp)
					{
						tmp = dblFuncValues[j];
					}
				}
				alpha[i] = tmp;
			}
			else if (payoffType == "avg")
			{
				tmp = 0.0;
				for (size_t j = 0; j < dim; j++)
				{
					tmp += dblFuncValues[j];
				}
				alpha[i] = tmp/static_cast<double>(dim);
			}
			else if (payoffType == "avgM")
			{
				tmp = 0.0;
				for (size_t j = 0; j < dim; j++)
				{
					tmp += dblFuncValues[j];
				}
				alpha[i] = max(((tmp/static_cast<double>(dim))-1.0), 0.0);
			}
			else
			{
				throw new application_exception("BlackScholesSolver::initGridWithEuroCallPayoff : An unknown payoff-type was specified!");
			}

			delete[] dblFuncValues;
		}

		OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::initGridWithEuroCallPayoff : A grid wasn't constructed before!");
	}
}

double BlackScholesSolver::getOptionPrice(std::vector<double>& evalPoint, DataVector& alpha)
{
	double result = 0.0;

	if (bGridConstructed)
	{
		OperationEval* myEval = myGrid->createOperationEval();
		result = myEval->eval(alpha, evalPoint);
		delete myEval;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::getOptionPrice : A grid wasn't constructed before!");
	}

	return result;
}

size_t BlackScholesSolver:: getNumberGridPoints()
{
	if (bGridConstructed)
	{
		return myGridStorage->size();
	}
	else
	{
		throw new application_exception("BlackScholesSolver::getNumberGridPoints : A grid wasn't constructed before!");
	}
}

size_t BlackScholesSolver:: getNumberDimensions()
{
	if (bGridConstructed)
	{
		return myGridStorage->dim();
	}
	else
	{
		throw new application_exception("BlackScholesSolver::getNumberDimensions : A grid wasn't constructed before!");
	}
}

double BlackScholesSolver::get1DEuroCallPayoffValue(double assetValue, double strike)
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

void BlackScholesSolver::solve1DAnalytic(std::vector< std::pair<double, double> >& premiums, double maxStock, double StockInc, double strike, double t)
{
	if (bStochasticDataAlloc)
	{
		double stock = 0.0;
		double vola = this->sigmas->get(0);
		StdNormalDistribution* myStdNDis = new StdNormalDistribution();

		for (stock = 0.0; stock <= maxStock; stock += StockInc)
		{
			double dOne = (log((stock/strike)) + ((this->r + (vola*vola*0.5))*(t)))/(vola*sqrt(t));
			double dTwo = dOne - (vola*sqrt(t));
			double prem = (stock*myStdNDis->getCumulativeDensity(dOne)) - (strike*myStdNDis->getCumulativeDensity(dTwo)*(exp((-1.0)*this->r*t)));

			premiums.push_back(std::make_pair(stock, prem));
		}

		delete myStdNDis;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::solve1DAnalytic : Stochastic parameters weren't set!");
	}
}

void BlackScholesSolver::print1DAnalytic(std::vector< std::pair<double, double> >& premiums, std::string tfilename)
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

void BlackScholesSolver::initScreen()
{
	myScreen = new ScreenOutput();
	myScreen->writeTitle("SGpp - Black Scholes Solver, 1.0.1", "Alexander Heinecke, (C) 2009-2010");
	myScreen->writeStartSolve("Multidimensional Black Scholes Solver");
}

void BlackScholesSolver::writeHelp()
{
	std::stringstream mySStream;

	if (myScreen != NULL)
	{
		mySStream << "Some instructions for the use of Black Scholes Solver:" << std::endl;
		mySStream << "------------------------------------------------------" << std::endl << std::endl;
		mySStream << "Available execution modes are:" << std::endl;
		mySStream << "	test1D		Solves a simple 1D example" << std::endl;
		mySStream << "	test2D		Solves a 2D example" << std::endl;
		mySStream << "	solveND		Solves a ND example" << std::endl;
		mySStream << "	solveBonn	Solves an option delivered in Bonn's format" << std::endl << std::endl;

		mySStream << "test1D" << std::endl << "------" << std::endl;
		mySStream << "the following options must be specified:" << std::endl;
		mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
		mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
		mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
		mySStream << "	strike: the strike of the call option" << std::endl;
		mySStream << "	r: the riskfree rate" << std::endl;
		mySStream << "	T: time to maturity" << std::endl;
		mySStream << "	dT: timestep size" << std::endl;
		mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
		mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
		mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
		mySStream << "	[no]animation: generate pictures for an animation" << std::endl;
		mySStream << std::endl;
		mySStream << "Example:" << std::endl;
		mySStream << "5 " << "bound.data stoch.data " << "65.0 " << "0.05 " << "1.0 " << "0.1 ImEul " << "400 " << "0.000001 " << "noanimation" << std::endl;
		mySStream << std::endl;
		mySStream << "Remark: This test generates following output files:" << std::endl;
		mySStream << "	payoff.gnuplot: the start condition" << std::endl;
		mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
		mySStream << "	analyticBS.gnuplot: the analytical solution" << std::endl;
		mySStream << std::endl << std::endl;

		mySStream << "test2D" << std::endl << "------" << std::endl;
		mySStream << "the following options must be specified:" << std::endl;
		mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
		mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
		mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
		mySStream << "	strike1: the strike 1 of the call option" << std::endl;
		mySStream << "	strike2: the strike 1 of the call option" << std::endl;
		mySStream << "	r: the riskfree rate" << std::endl;
		mySStream << "	T: time to maturity" << std::endl;
		mySStream << "	dT: timestep size" << std::endl;
		mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
		mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
		mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
		mySStream << "	[no]animation: generate pictures for an animation" << std::endl;
		mySStream << std::endl;
		mySStream << "Example:" << std::endl;
		mySStream << "5 " << "bound.data stoch.data " << "65.0 55.0 " << "0.05 " << "1.0 " << "0.1 ImEul " << "400 " << "0.000001 " << "noanimation" << std::endl;
		mySStream << std::endl;
		mySStream << "Remark: This test generates following output files:" << std::endl;
		mySStream << "	payoff.gnuplot: the start condition" << std::endl;
		mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
		mySStream << std::endl << std::endl;

		mySStream << "solveND" << std::endl << "------" << std::endl;
		mySStream << "the following options must be specified:" << std::endl;
		mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
		mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
		mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
		mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
		mySStream << "	file_Strikes: file containing strikes of Europ. Call" << std::endl;
		mySStream << "	payoff_func: function for n-d payoff: max, avg/avgM" << std::endl;
		mySStream << "	r: the riskfree rate" << std::endl;
		mySStream << "	T: time to maturity" << std::endl;
		mySStream << "	dT: timestep size" << std::endl;
		mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
		mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
		mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
		mySStream << std::endl;
		mySStream << "Example:" << std::endl;
		mySStream << "3 5 " << "bound.data stoch.data strike.data max "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 " << std::endl;
		mySStream << std::endl;
		mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
		mySStream << "	payoff.gnuplot: the start condition" << std::endl;
		mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
		mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
		mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
		mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
		mySStream << std::endl << std::endl;

		mySStream << "solveNDadapt" << std::endl << "------" << std::endl;
		mySStream << "the following options must be specified:" << std::endl;
		mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
		mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
		mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
		mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
		mySStream << "	file_Strikes: file containing strikes of Europ. Call" << std::endl;
		mySStream << "	payoff_func: function for n-d payoff: max, avg/avgM" << std::endl;
		mySStream << "	r: the riskfree rate" << std::endl;
		mySStream << "	T: time to maturity" << std::endl;
		mySStream << "	dT: timestep size" << std::endl;
		mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
		mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
		mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
		mySStream << "	Adapt-Initial-Refinement: Number of Initial" << std::endl;
		mySStream << "			Refinements" << std::endl;
		mySStream << "	Adapt-Initial-Distance: determines the distance" << std::endl;
		mySStream << "			a grid point must have from @money to" << std::endl;
		mySStream << "			by refined" << std::endl;
		mySStream << std::endl;
		mySStream << "Example:" << std::endl;
		mySStream << "3 5 " << "bound.data stoch.data strike.data max "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 5 0.5" << std::endl;
		mySStream << std::endl;
		mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
		mySStream << "	payoff.gnuplot: the start condition" << std::endl;
		mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
		mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
		mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
		mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
		mySStream << std::endl << std::endl;

		mySStream << "solveBonn" << std::endl << "---------" << std::endl;
		mySStream << "the following options must be specified:" << std::endl;
		mySStream << "	file_grid_in: file the specifies the unsolved grid" << std::endl;
		mySStream << "	file_grid_out: file that contains the solved grid when finished" << std::endl;
		mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
		mySStream << "	r: the riskfree rate" << std::endl;
		mySStream << "	T: time to maturity" << std::endl;
		mySStream << "	dT: timestep size" << std::endl;
		mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
		mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
		mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
		mySStream << std::endl;
		mySStream << "Example:" << std::endl;
		mySStream << "grid.in grid.out " << "stoch.data " << "0.05 " << "1.0 " << "0.1 ImEul " << "400 " << "0.000001 " << std::endl;

		mySStream << std::endl << std::endl;
		myScreen->writeHelp(mySStream.str());
	}
}

}
