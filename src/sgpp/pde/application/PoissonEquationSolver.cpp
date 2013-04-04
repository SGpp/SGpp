/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/application/PoissonEquationSolver.hpp"
#include "pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "base/grid/Grid.hpp"
#include "base/exception/application_exception.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "stdlib.h"
#include <sstream>
#include <fstream>
using namespace sg::solver;
using namespace sg::base;

namespace sg
{
namespace pde
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

void PoissonEquationSolver::constructGrid(BoundingBox& BoundingBox, int level)
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

		OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::initGridWithSmoothHeat : A grid wasn't constructed before!");
	}
}

void PoissonEquationSolver::initGridWithSmoothHeatFullDomain(DataVector& alpha, double mu, double sigma, double factor)
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

		OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("HeatEquationSolver::initGridWithSmoothHeatFullDomain : A grid wasn't constructed before!");
	}
}

void PoissonEquationSolver::initGridWithExpHeat(DataVector& alpha, double factor)
{
	if (this->bGridConstructed)
	{
		double tmp;
		double* dblFuncValues = new double[this->dim];
		double* rightBound = new double[this->dim];

		BoundingBox* tmpBB = this->myGrid->getBoundingBox();

		for (size_t j = 0; j < this->dim; j++)
		{
			rightBound[j] = (tmpBB->getBoundary(j)).rightBoundary;
		}

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
			std::stringstream coordsStream(coords);
			bool isInner = true;
			tmp = 0.0;

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
					tmp *= exp((dblFuncValues[j]-rightBound[j])*factor);
				}
			}
			else
			{
				tmp = 0.0;
			}

			alpha[i] = tmp;
		}

		delete[] dblFuncValues;

		OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("PoissonEquationSolver::initGridWithExpHeat : A grid wasn't constructed before!");
	}
}

void PoissonEquationSolver::initGridWithExpHeatFullDomain(DataVector& alpha, double factor)
{
	if (this->bGridConstructed)
	{
		double tmp;
		double* dblFuncValues = new double[this->dim];
		double* rightBound = new double[this->dim];

		BoundingBox* tmpBB = this->myGrid->getBoundingBox();

		for (size_t j = 0; j < this->dim; j++)
		{
			rightBound[j] = (tmpBB->getBoundary(j)).rightBoundary;
		}

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
			std::stringstream coordsStream(coords);
			tmp = 0.0;

			for (size_t j = 0; j < this->dim; j++)
			{
				coordsStream >> tmp;

				dblFuncValues[j] = tmp;
			}

			tmp = 1.0;
			for (size_t j = 0; j < this->dim; j++)
			{
				tmp *= exp((dblFuncValues[j]-rightBound[j])*factor);
			}

			alpha[i] = tmp;
		}

		delete[] dblFuncValues;

		OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
		myHierarchisation->doHierarchisation(alpha);
		delete myHierarchisation;
	}
	else
	{
		throw new application_exception("PoissonEquationSolver::initGridWithExpHeat : A grid wasn't constructed before!");
	}
}

void PoissonEquationSolver::storeInnerMatrix(std::string tFilename)
{
	DataVector rhs(this->myGrid->getSize());
	rhs.setAll(0.0);
	PoissonEquationEllipticPDESolverSystemDirichlet* mySystem = new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), rhs);
	SGppStopwatch* myStopwatch = new SGppStopwatch();

	std::string mtx = "";

	std::cout << "Generating matrix in MatrixMarket format..." << std::endl;
	myStopwatch->start();
	mySystem->getInnerMatrix(mtx);

	std::ofstream outfile(tFilename.c_str());
	outfile << mtx;
	outfile.close();
	std::cout << "Generating matrix in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;

	delete myStopwatch;
	delete mySystem;
}

void PoissonEquationSolver::storeInnerMatrixDiagonal(std::string tFilename)
{
	DataVector rhs(this->myGrid->getSize());
	rhs.setAll(0.0);
	PoissonEquationEllipticPDESolverSystemDirichlet* mySystem = new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), rhs);
	SGppStopwatch* myStopwatch = new SGppStopwatch();

	std::string mtx = "";

	std::cout << "Generating diagonal matrix in MatrixMarket format..." << std::endl;
	myStopwatch->start();
	mySystem->getInnerMatrixDiagonal(mtx);

	std::ofstream outfile(tFilename.c_str());
	outfile << mtx;
	outfile.close();
	std::cout << "Generating diagonal matrix in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;

	delete myStopwatch;
	delete mySystem;
}

void PoissonEquationSolver::storeInnerMatrixDiagonalRowSum(std::string tFilename)
{
	DataVector rhs(this->myGrid->getSize());
	rhs.setAll(0.0);
	PoissonEquationEllipticPDESolverSystemDirichlet* mySystem = new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), rhs);
	SGppStopwatch* myStopwatch = new SGppStopwatch();

	std::string mtx = "";

	std::cout << "Generating row sum diagonal matrix in MatrixMarket format..." << std::endl;
	myStopwatch->start();
	mySystem->getInnerMatrixDiagonalRowSum(mtx);

	std::ofstream outfile(tFilename.c_str());
	outfile << mtx;
	outfile.close();
	std::cout << "Generating row sum diagonal matrix in MatrixMarket format... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;

	delete myStopwatch;
	delete mySystem;
}

void PoissonEquationSolver::storeInnerRHS(DataVector& alpha, std::string tFilename)
{
	SGppStopwatch* myStopwatch = new SGppStopwatch();
	PoissonEquationEllipticPDESolverSystemDirichlet* mySystem = new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), alpha);

	std::cout << "Exporting inner right-hand-side..." << std::endl;
	myStopwatch->start();
	DataVector* rhs_inner = mySystem->generateRHS();

	size_t nCoefs = rhs_inner->getSize();
	std::ofstream outfile(tFilename.c_str());
	for (size_t i = 0; i < nCoefs; i++)
	{
		outfile << std::scientific << rhs_inner->get(i) << std::endl;
	}
	outfile.close();
	std::cout << "Exporting inner right-hand-side... DONE! (" << myStopwatch->stop() << " s)" << std::endl << std::endl << std::endl;

	delete mySystem; // rhs_inner are allocated and freed here!!
	delete myStopwatch;
}

void PoissonEquationSolver::storeInnerSolution(DataVector& alpha, size_t maxCGIterations, double epsilonCG, std::string tFilename)
{
	ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
	PoissonEquationEllipticPDESolverSystemDirichlet* mySystem = new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), alpha);

	std::cout << "Exporting inner solution..." << std::endl;

	DataVector* alpha_solve = mySystem->getGridCoefficientsForCG();
	DataVector* rhs_solve = mySystem->generateRHS();

	myCG->solve(*mySystem, *alpha_solve, *rhs_solve, true, false, 0.0);

	size_t nCoefs = alpha_solve->getSize();
	std::ofstream outfile(tFilename.c_str());
	for (size_t i = 0; i < nCoefs; i++)
	{
		outfile << std::scientific << alpha_solve->get(i) << std::endl;
	}
	outfile.close();

	std::cout << "Exporting inner solution... DONE!" << std::endl;

	delete myCG;
	delete mySystem; // alpha_solver and rhs_solve are allocated and freed here!!
}

void PoissonEquationSolver::initScreen()
{
	this->myScreen = new ScreenOutput();
	this->myScreen->writeTitle("SGpp - Poisson Equation Solver, 1.0.0", "Alexander Heinecke, (C) 2009-2011");
}

}
}
