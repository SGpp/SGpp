/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#include "algorithm/pde/BlackScholesTimestepMatrix.hpp"
#include "application/pde/BlackScholesSolver.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "grid/Grid.hpp"
#include "stdlib.h"
#include <sstream>

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

void BlackScholesSolver::constructGrid(std::string tfilename, DataVector& emptyAlpha, bool& ishierarchized)
{
	IOToolBonnSG* myImporter = new IOToolBonnSG();
	std::string serGrid;

	// test if emptyAlpha is really empty
	if (emptyAlpha.getSize() != 0)
	{
		// @todo (heinecke) thrown an application exception
		return;
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
		// @todo (heinecke) through an application exception here
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
		Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, maxCGIterations, epsilonCG, generateAnimation, numEvalsAnimation, myScreen);
		BlackScholesTimestepMatrix* myBSMatrix = new BlackScholesTimestepMatrix(*myGrid, *this->mus, *this->sigmas, *this->rhos, r, timestepsize, "ExEul");

		myEuler->solve(*myBSMatrix, alpha, verbose);

		delete myBSMatrix;
		delete myEuler;
	}
	else
	{
		// @todo (heinecke) through an application exception here
	}
}

void BlackScholesSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (bGridConstructed && bStochasticDataAlloc)
	{
		Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, maxCGIterations, epsilonCG, generateAnimation, numEvalsAnimation, myScreen);
		BlackScholesTimestepMatrix* myBSMatrix = new BlackScholesTimestepMatrix(*myGrid, *this->mus, *this->sigmas, *this->rhos, r, timestepsize, "ImEul");

		myEuler->solve(*myBSMatrix, alpha, verbose);

		delete myBSMatrix;
		delete myEuler;
	}
	else
	{
		// @todo (heinecke) through an application exception here
	}
}

void BlackScholesSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha)
{
	if (bGridConstructed && bStochasticDataAlloc)
	{
		CrankNicolson* myCN = new CrankNicolson(numTimesteps, timestepsize, maxCGIterations, epsilonCG);
		BlackScholesTimestepMatrix* myBSMatrix = new BlackScholesTimestepMatrix(*myGrid, *this->mus, *this->sigmas, *this->rhos, r, timestepsize, "CrNic");

		myCN->solve(*myBSMatrix, alpha, true);

		delete myBSMatrix;
		delete myCN;
	}
	else
	{
		// @todo (heinecke) through an application exception here
	}
}

void BlackScholesSolver::printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename)
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void BlackScholesSolver::initGridWithPayoff(DataVector& alpha, double* strike)
{
	double tmp;
	double tmp2;

	if (bGridConstructed)
	{
		if (dim == 1)
		{
			for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
			{
				tmp = atof(myGridStorage->get(i)->getCoordsStringBB(*myBoundingBox).c_str());
				alpha[i] = get1DPayoffValue(tmp, strike[0]);
			}

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else if (dim == 2)
		{
			for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
			{
					std::string coords = myGridStorage->get(i)->getCoordsStringBB(*myBoundingBox);
					std::stringstream coordsStream(coords);

					coordsStream >> tmp;
					coordsStream >> tmp2;

					alpha[i] = max(get1DPayoffValue(tmp, strike[0]),get1DPayoffValue(tmp2, strike[1]));
			}

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{
			// @todo (heinecke) thrown an application exception
		}
	}
	else
	{
		// @todo (heinecke) throw an application exception
	}
}

size_t BlackScholesSolver:: getNumberGridPoints()
{
	if (bGridConstructed)
	{
		return myGridStorage->size();
	}
	else
	{
		// @todo (heinecke) throw an application exception
		return 0;
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
		// @todo (heinecke) throw an application exception
		return 0;
	}
}

double BlackScholesSolver::get1DPayoffValue(double assetValue, double strike)
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
		// @todo (heinecke) throw an application exception
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
	myScreen->writeTitle("SGpp - Black Scholes Solver, 1.0.0 alpha", "Alexander Heinecke, (C) 2009");
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
		mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
		mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
		mySStream << "	[no]animation: generate pictures for an animation" << std::endl;
		mySStream << std::endl;
		mySStream << "Example:" << std::endl;
		mySStream << "5 " << "stoch.data " << "bound.data " << "65.0 " << "0.05 " << "1.0 " << "0.1 " << "400 " << "0.000001 " << "noanimation" << std::endl;
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
		mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
		mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
		mySStream << "	[no]animation: generate pictures for an animation" << std::endl;
		mySStream << std::endl;
		mySStream << "Example:" << std::endl;
		mySStream << "5 " << "stoch.data " << "bound.data " << "65.0 55.0 " << "0.05 " << "1.0 " << "0.1 " << "400 " << "0.000001 " << "noanimation" << std::endl;
		mySStream << std::endl;
		mySStream << "Remark: This test generates following output files:" << std::endl;
		mySStream << "	payoff.gnuplot: the start condition" << std::endl;
		mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
		mySStream << std::endl << std::endl;

		mySStream << "solveBonn" << std::endl << "---------" << std::endl;
		mySStream << "the following options must be specified:" << std::endl;
		mySStream << "	file_grid_in: file the specifies the unsolved grid" << std::endl;
		mySStream << "	file_grid_out: file that contains the solved grid when finished" << std::endl;
		mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
		mySStream << "	r: the riskfree rate" << std::endl;
		mySStream << "	T: time to maturity" << std::endl;
		mySStream << "	dT: timestep size" << std::endl;
		mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
		mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
		mySStream << std::endl;
		mySStream << "Example:" << std::endl;
		mySStream << "grid.in grid.out " << "stoch.data " << "0.05 " << "1.0 " << "0.1 " << "400 " << "0.000001 " << std::endl;

		mySStream << std::endl << std::endl;
		myScreen->writeHelp(mySStream.str());
	}
}

}
