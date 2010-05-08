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

void BlackScholesSolver::deleteGrid()
{
	if (bGridConstructed)
	{
		delete myGrid;
		bGridConstructed = false;
		myBoundingBox = NULL;
		myGridStorage = NULL;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::deleteGrid : The grid wasn't initialized before!");
	}
}

void BlackScholesSolver::setGrid(std::string serializedGrid)
{
	if (bGridConstructed)
	{
		delete myGrid;
		bGridConstructed = false;
		myBoundingBox = NULL;
		myGridStorage = NULL;
	}

	myGrid = Grid::unserialize(serializedGrid);

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	bGridConstructed = true;
}

std::string BlackScholesSolver::getGrid()
{
	std::string gridSer = "";

	if (bGridConstructed)
	{
		// Serialize the grid
		myGrid->serialize(gridSer);
	}
	else
	{
		throw new application_exception("BlackScholesSolver::getGrid : The grid wasn't initialized before!");
	}

	return gridSer;
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

void BlackScholesSolver::getCuboidEvalPoints(std::vector<DataVector>& evalPoints, DataVector& curPoint, std::vector<double>& center, double size, size_t points, size_t curDim)
{
	if (curDim == 0)
	{
		if (points > 1)
		{
			for (size_t i = 0; i < points; i++)
			{
				curPoint.set(curDim, min(max(center[curDim]-(myBoundingBox->getIntervalWidth(curDim)*size)+
						((myBoundingBox->getIntervalWidth(curDim)*size*2/(points-1))*static_cast<double>(i)),
						myBoundingBox->getIntervalOffset(curDim)),
						myBoundingBox->getIntervalOffset(curDim)+myBoundingBox->getIntervalWidth(curDim)));

				evalPoints.push_back(curPoint);
			}
		}
		else
		{
			curPoint.set(curDim, center[curDim]);

			evalPoints.push_back(curPoint);
		}
	}
	else
	{
		if (points > 1)
		{
			for (size_t i = 0; i < points; i++)
			{
				curPoint.set(curDim, min(max(center[curDim]-(myBoundingBox->getIntervalWidth(curDim)*size)+
						((myBoundingBox->getIntervalWidth(curDim)*size*2/(points-1))*static_cast<double>(i)),
						myBoundingBox->getIntervalOffset(curDim)),
						myBoundingBox->getIntervalOffset(curDim)+myBoundingBox->getIntervalWidth(curDim)));

				getCuboidEvalPoints(evalPoints, curPoint, center, size, points, curDim-1);
			}
		}
		else
		{
			curPoint.set(curDim, center[curDim]);

			getCuboidEvalPoints(evalPoints, curPoint, center, size, points, curDim-1);
		}
	}
}

void BlackScholesSolver::getOptionPricesCuboid(DataVector& alpha, DataVector& OptionPrices, DataVector& EvaluationPoints)
{
	if (bGridConstructed)
	{
		if (OptionPrices.getSize() != EvaluationPoints.getSize())
		{
			throw new application_exception("BlackScholesSolver::getOptionPricesCuboid : The size of the price vector doesn't match the size of the evaluation points' vector!");
		}

		OperationB* myOpB = myGrid->createOperationB();
		myOpB->multTranspose(alpha, EvaluationPoints, OptionPrices);
		delete myOpB;
	}
	else
	{
		throw new application_exception("BlackScholesSolver::getOptionPricesCuboid : A grid wasn't constructed before!");
	}
}

void BlackScholesSolver::getEvaluationCuboid(DataVector& EvaluationPoints, std::vector<double>& center, double size, size_t points)
{
	std::vector<DataVector> evalPoints;
	DataVector curPoint(getNumberDimensions());

	getCuboidEvalPoints(evalPoints, curPoint, center, size, points, getNumberDimensions()-1);

	size_t numEvalPoints = evalPoints.size();
	EvaluationPoints.resize(numEvalPoints);

	for (size_t i = 0; i < numEvalPoints; i++)
	{
		EvaluationPoints.setRow(i, evalPoints[i]);
	}
}

size_t BlackScholesSolver::getNumberGridPoints()
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

size_t BlackScholesSolver::getNumberInnerGridPoints()
{
	if (bGridConstructed)
	{
		return myGridStorage->getNumInnerPoints();
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
	myScreen->writeTitle("SGpp - Black Scholes Solver, 1.1.0", "Alexander Heinecke, (C) 2009-2010");
	myScreen->writeStartSolve("Multidimensional Black Scholes Solver");
}

}
