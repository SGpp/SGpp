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

#include "algorithm/pde/HeatEquationODESolverSystem.hpp"
#include "application/pde/HeatEquationSolver.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include "stdlib.h"
#include <sstream>

namespace sg
{

HeatEquationSolver::HeatEquationSolver()
{
	bGridConstructed = false;
	myScreen = NULL;
}

HeatEquationSolver::~HeatEquationSolver()
{
	if (bGridConstructed)
	{
		delete myGrid;
	}
	if (myScreen != NULL)
	{
		delete myScreen;
	}
}

void HeatEquationSolver::constructGrid(BoundingBox& BoundingBox, size_t level, bool useBoundary)
{
	dim = BoundingBox.getDimensions();
	levels = level;

	if (useBoundary)
	{
		myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);
	}
	else
	{
		myGrid = new LinearGrid(dim);
	}

	GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(levels);
	delete myGenerator;

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	bGridConstructed = true;
}

void HeatEquationSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (bGridConstructed)
	{
		myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*myGrid, alpha, a, timestepsize, "ExEul");

		myEuler->solve(*myCG, *myHESolver, verbose);

		delete myHESolver;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::solveExplicitEuler : A grid wasn't constructed before!");
	}
}

void HeatEquationSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (bGridConstructed)
	{
		myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*myGrid, alpha, a, timestepsize, "ImEul");

		myEuler->solve(*myCG, *myHESolver, verbose);

		delete myHESolver;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
	}
}

void HeatEquationSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha)
{
	if (bGridConstructed)
	{
		myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		CrankNicolson* myCN = new CrankNicolson(numTimesteps, timestepsize);
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*myGrid, alpha, a, timestepsize, "CrNic");

		myCN->solve(*myCG, *myHESolver, false);

		delete myHESolver;
		delete myCG;
		delete myCN;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::solveCrankNicolson : A grid wasn't constructed before!");
	}
}

void HeatEquationSolver::printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename)
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void HeatEquationSolver::initGridWithSingleHeat(DataVector& alpha, double heat)
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

				if (tmp == 0.5)
				{
					alpha[i] = heat;
				}
				else
				{
					alpha[i] = 0.0;
				}
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

					if (tmp == 0.5 && tmp2 == 0.5)
					{
						alpha[i] = heat;
					}
					else
					{
						alpha[i] = 0.0;
					}
			}

			//std::cout << alpha.toString() << std::endl;

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);

			//std::cout << alpha.toString() << std::endl;

			//myHierarchisation->doDehierarchisation(alpha);

			//std::cout << alpha.toString() << std::endl;
			delete myHierarchisation;
		}
		else
		{
			throw new application_exception("HeatEquationSolver::initGridWithSingleHeat : The constructed grid has more than two dimensions!");
		}
	}
	else
	{
		throw new application_exception("HeatEquationSolver::initGridWithSingleHeat : A grid wasn't constructed before!");
	}
}

size_t HeatEquationSolver:: getNumberGridPoints()
{
	if (bGridConstructed)
	{
		return myGridStorage->size();
	}
	else
	{
		throw new application_exception("HeatEquationSolver::getNumberGridPoints : A grid wasn't constructed before!");
	}
}

void HeatEquationSolver::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor)
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

				alpha[i] = factor*(1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp-mu)/sigma)*((tmp-mu)/sigma));
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

				alpha[i] = factor*factor*((1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp-mu)/sigma)*((tmp-mu)/sigma))) * ((1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp2-mu)/sigma)*((tmp2-mu)/sigma)));
			}

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{
			throw new application_exception("HeatEquationSolver::initGridWithSmoothHeat : The constructed grid has more than two dimensions!");
		}
	}
	else
	{
		throw new application_exception("HeatEquationSolver::initGridWithSmoothHeat : A grid wasn't constructed before!");
	}
}

void HeatEquationSolver::initGridWithConstantHeat(DataVector& alpha, double constHeat)
{
	double tmp;
	//double tmp2;

	if (bGridConstructed)
	{
		if (dim == 1)
		{
			for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
			{
				tmp = atof(myGridStorage->get(i)->getCoordsStringBB(*myBoundingBox).c_str());

				alpha[i] = constHeat;
			}

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{
			throw new application_exception("HeatEquationSolver::initGridWithConstantHeat : The constructed grid has more than one dimension!");
		}
	}
	else
	{
		throw new application_exception("HeatEquationSolver::initGridWithConstantHeat : A grid wasn't constructed before!");
	}
}

void HeatEquationSolver::initScreen()
{
	myScreen = new ScreenOutput();
	myScreen->writeTitle("SGpp - Heat Equation Solver, 1.0.0", "Alexander Heinecke, (C) 2009-2010");
}

}
