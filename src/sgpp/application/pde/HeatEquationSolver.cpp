/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

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

HeatEquationSolver::HeatEquationSolver() : ParabolicPDESolver()
{
	this->bGridConstructed = false;
	this->myScreen = NULL;
}

HeatEquationSolver::~HeatEquationSolver()
{
	if (this->myScreen != NULL)
	{
		delete this->myScreen;
	}
}

void HeatEquationSolver::constructGrid(BoundingBox& BoundingBox, size_t level)
{
	this->dim = BoundingBox.getDimensions();
	this->levels = level;

	this->myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

	GridGenerator* myGenerator = this->myGrid->createGridGenerator();
	myGenerator->regular(this->levels);
	delete myGenerator;

	this->myBoundingBox = this->myGrid->getBoundingBox();
	this->myGridStorage = this->myGrid->getStorage();

	this->bGridConstructed = true;
}

void HeatEquationSolver::setHeatCoefficient(double a)
{
	this->a = a;
}

void HeatEquationSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed)
	{
		this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ExEul");

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

void HeatEquationSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed)
	{
		this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "ImEul");

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

void HeatEquationSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul)
{
	if (this->bGridConstructed)
	{
		this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
		HeatEquationODESolverSystem* myHESolver = new HeatEquationODESolverSystem(*this->myGrid, alpha, this->a, timestepsize, "CrNic");

		size_t numCNSteps;
		size_t numIESteps;

		numCNSteps = numTimesteps;
		if (numTimesteps > NumImEul)
		{
			numCNSteps = numTimesteps - NumImEul;
		}
		numIESteps = NumImEul;

		Euler* myEuler = new Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
		CrankNicolson* myCN = new CrankNicolson(numCNSteps, timestepsize);

		if (numIESteps > 0)
		{
			myEuler->solve(*myCG, *myHESolver, false);
		}

		myCN->solve(*myCG, *myHESolver, false);

		delete myHESolver;
		delete myCG;
		delete myCN;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::solveCrankNicolson : A grid wasn't constructed before!");
	}
}

void HeatEquationSolver::initGridWithSingleHeat(DataVector& alpha, double heat)
{
	double tmp;
	double tmp2;

	if (this->bGridConstructed)
	{
		if (this->dim == 1)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				tmp = atof(this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox).c_str());

				if (tmp == 0.5)
				{
					alpha[i] = heat;
				}
				else
				{
					alpha[i] = 0.0;
				}
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else if (this->dim == 2)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
					std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
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

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
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

void HeatEquationSolver::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor)
{
	double tmp;
	double tmp2;

	if (this->bGridConstructed)
	{
		if (this->dim == 1)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				tmp = atof(this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox).c_str());

				alpha[i] = factor*(1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp-mu)/sigma)*((tmp-mu)/sigma));
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else if (this->dim == 2)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
				std::stringstream coordsStream(coords);

				coordsStream >> tmp;
				coordsStream >> tmp2;

				alpha[i] = factor*factor*((1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp-mu)/sigma)*((tmp-mu)/sigma))) * ((1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp2-mu)/sigma)*((tmp2-mu)/sigma)));
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
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

	if (this->bGridConstructed)
	{
		if (this->dim == 1)
		{
			for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
			{
				tmp = atof(this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox).c_str());

				alpha[i] = constHeat;
			}

			OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
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
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Heat Equation Solver, 1.0.0", "Alexander Heinecke, (C) 2009-2010");
}

}
