/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/MPI/SGppMPITools.hpp"
#include "solver/sle/ConjugateGradientsMPI.hpp"
#include "application/pde/HeatEquationSolverMPI.hpp"
#include "algorithm/pde/HeatEquationParabolicPDESolverSystemParallelMPI.hpp"

#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include "stdlib.h"
#include <sstream>

namespace sg
{

HeatEquationSolverMPI::HeatEquationSolverMPI() : ParabolicPDESolver()
{
	this->bGridConstructed = false;
	this->myScreen = NULL;
}

HeatEquationSolverMPI::~HeatEquationSolverMPI()
{
	if (this->myScreen != NULL)
	{
		delete this->myScreen;
	}
}

void HeatEquationSolverMPI::constructGrid(BoundingBox& BoundingBox, size_t level)
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

void HeatEquationSolverMPI::setHeatCoefficient(double a)
{
	this->a = a;
}

void HeatEquationSolverMPI::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed)
	{
		if (this->myScreen != NULL)
		{
			this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		}
		double dNeededTime;
		Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
		ConjugateGradientsMPI* myCG = new ConjugateGradientsMPI(maxCGIterations, epsilonCG);
		HeatEquationParabolicPDESolverSystemParallelMPI* myHESolver = new HeatEquationParabolicPDESolverSystemParallelMPI(*this->myGrid, alpha, this->a, timestepsize, "ExEul");
		SGppStopwatch* myStopwatch = new SGppStopwatch();

		myStopwatch->start();
		myEuler->solve(*myCG, *myHESolver, verbose);
		dNeededTime = myStopwatch->stop();

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myStopwatch;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::solveExplicitEuler : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverMPI::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation)
{
	if (this->bGridConstructed)
	{
		if (this->myScreen != NULL)
		{
			this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		}

		double dNeededTime;
		Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
		ConjugateGradientsMPI* myCG = new ConjugateGradientsMPI(maxCGIterations, epsilonCG);
		HeatEquationParabolicPDESolverSystemParallelMPI* myHESolver = new HeatEquationParabolicPDESolverSystemParallelMPI(*this->myGrid, alpha, this->a, timestepsize, "ImEul");
		SGppStopwatch* myStopwatch = new SGppStopwatch();

		myStopwatch->start();
		myEuler->solve(*myCG, *myHESolver, verbose);
		dNeededTime = myStopwatch->stop();

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myStopwatch;
		delete myHESolver;
		delete myCG;
		delete myEuler;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverMPI::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul)
{
	if (this->bGridConstructed)
	{
		if (this->myScreen != NULL)
		{
			this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
		}
		double dNeededTime;
		ConjugateGradientsMPI* myCG = new ConjugateGradientsMPI(maxCGIterations, epsilonCG);
		HeatEquationParabolicPDESolverSystemParallelMPI* myHESolver = new HeatEquationParabolicPDESolverSystemParallelMPI(*this->myGrid, alpha, this->a, timestepsize, "CrNic");
		SGppStopwatch* myStopwatch = new SGppStopwatch();

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

		myStopwatch->start();
		if (numIESteps > 0)
		{
			myEuler->solve(*myCG, *myHESolver, false);
		}

		myCN->solve(*myCG, *myHESolver, false);
		dNeededTime = myStopwatch->stop();

		if (this->myScreen != NULL)
		{
			std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
			this->myScreen->writeEmptyLines(2);
		}

		delete myStopwatch;
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

void HeatEquationSolverMPI::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor)
{
	if (this->bGridConstructed)
	{
		double tmp;
		double* dblFuncValues = new double[this->dim];

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
			std::stringstream coordsStream(coords);

			for (size_t j = 0; j < this->dim; j++)
			{
				coordsStream >> tmp;

				dblFuncValues[j] = tmp;
			}

			tmp = 1.0;
			for (size_t j = 0; j < this->dim; j++)
			{
				tmp *=  factor*factor*((1.0/(sigma*2.0*3.145))*exp((-0.5)*((dblFuncValues[j]-mu)/sigma)*((dblFuncValues[j]-mu)/sigma)));
			}

			alpha[i] = tmp;
		}

		delete[] dblFuncValues;

		OperationHierarchisation* myHierarchisation = this->myGrid->createOperationHierarchisation();
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::initGridWithSmoothHeat : A grid wasn't constructed before!");
	}
}

void HeatEquationSolverMPI::initScreen()
{
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Heat Equation Solver, 1.0.0", "Alexander Heinecke, (C) 2009-2011");
}

}
