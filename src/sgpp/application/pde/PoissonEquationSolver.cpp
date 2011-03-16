/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "application/pde/PoissonEquationSolver.hpp"
#include "algorithm/pde/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include "tools/common/SGppStopwatch.hpp"
#include "stdlib.h"
#include <sstream>
using namespace sg::base;

namespace sg
{

PoissonEquationSolver::PoissonEquationSolver() : EllipticPDESolver()
{
	this->bGridConstructed = false;
	this->myScreen = NULL;
}

PoissonEquationSolver::~PoissonEquationSolver()
{
	if (this->myScreen != NULL)
	{
		delete this->myScreen;
	}
}

void PoissonEquationSolver::constructGrid(BoundingBox& BoundingBox, size_t level)
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

void PoissonEquationSolver::solvePDE(DataVector& alpha, DataVector& rhs, size_t maxCGIterations, double epsilonCG, bool verbose)
{
	double dTimeAlpha = 0.0;
	double dTimeRHS = 0.0;
	double dTimeSolver = 0.0;

	SGppStopwatch* myStopwatch = new SGppStopwatch();
	ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
	PoissonEquationEllipticPDESolverSystemDirichlet* mySystem = new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), rhs);

	std::cout << "Gridpoints (complete grid): " << mySystem->getNumGridPointsComplete() << std::endl;
	std::cout << "Gridpoints (inner grid): " << mySystem->getNumGridPointsInner() << std::endl << std::endl << std::endl;

	myStopwatch->start();
	DataVector* alpha_solve = mySystem->getGridCoefficientsForCG();
	dTimeAlpha = myStopwatch->stop();
	std::cout << "coefficients has been initialized for solving!" << std::endl;
	myStopwatch->start();
	DataVector* rhs_solve = mySystem->generateRHS();
	dTimeRHS = myStopwatch->stop();
	std::cout << "right hand side has been initialized for solving!" << std::endl << std::endl << std::endl;

	myStopwatch->start();
	myCG->solve(*mySystem, *alpha_solve, *rhs_solve, true, verbose, 0.0);

	// Copy result into coefficient vector of the boundary grid
	mySystem->getSolutionBoundGrid(alpha, *alpha_solve);
	dTimeSolver = myStopwatch->stop();

	std::cout << std::endl << std::endl;
	std::cout << "Gridpoints (complete grid): " << mySystem->getNumGridPointsComplete() << std::endl;
	std::cout << "Gridpoints (inner grid): " << mySystem->getNumGridPointsInner() << std::endl << std::endl << std::endl;

	std::cout << "Timings for solving Poisson Equation" << std::endl;
	std::cout << "------------------------------------" << std::endl;
	std::cout << "Time for creating CG coeffs: " << dTimeAlpha << std::endl;
	std::cout << "Time for creating RHS: " << dTimeRHS << std::endl;
	std::cout << "Time for solving: " << dTimeSolver << std::endl << std::endl;
	std::cout << "Time: " << dTimeAlpha + dTimeRHS + dTimeSolver << std::endl << std::endl << std::endl;

	delete myCG;
	delete mySystem; // alpha_solver and rhs_solve are allocated and freed here!!
	delete myStopwatch;
}

void PoissonEquationSolver::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor)
{
	if (this->bGridConstructed)
	{
		double tmp;
		double* dblFuncValues = new double[this->dim];

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
			std::stringstream coordsStream(coords);
			bool isInner = true;

			for (size_t j = 0; j < this->dim; j++)
			{
				coordsStream >> tmp;

				// determine if a grid point is an inner grid point
				if ((tmp != this->myBoundingBox->getBoundary(j).leftBoundary && tmp != this->myBoundingBox->getBoundary(j).rightBoundary))
				{
					// Nothtin to do, test is that qay hence == for floating point values is unsave
				}
				else
				{
					isInner = false;
				}

				dblFuncValues[j] = tmp;
			}

			if (isInner == false)
			{
				tmp = 1.0;
				for (size_t j = 0; j < this->dim; j++)
				{
					tmp *=  factor*factor*((1.0/(sigma*2.0*3.145))*exp((-0.5)*((dblFuncValues[j]-mu)/sigma)*((dblFuncValues[j]-mu)/sigma)));
				}
			}
			else
			{
				tmp = 0.0;
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

void PoissonEquationSolver::initScreen()
{
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Poisson Equation Solver, 1.0.0", "Alexander Heinecke, (C) 2009-2011");
}

}
