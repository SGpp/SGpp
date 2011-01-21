/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "application/pde/PoissonEquationSolver.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include "stdlib.h"
#include <sstream>

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
	ConjugateGradients* myCG = new ConjugateGradients(maxCGIterations, epsilonCG);
	DirichletGridConverter* myGridConverter =  new DirichletGridConverter();
	DirichletUpdateVector* myBoundaryUpdate = new DirichletUpdateVector(this->myGrid->getStorage());
	OperationMatrix* L_complete = this->myGrid->createOperationLaplace();
	Grid* InnerGrid;
	DataVector* alpha_inner;
	DataVector* rhs_inner;
	DataVector alpha_complete(rhs);
	DataVector tmp_complete(rhs);

	myBoundaryUpdate->setInnerPointsToZero(alpha_complete);
	L_complete->mult(alpha_complete, tmp_complete);

	myGridConverter->buildInnerGridWithCoefs(*this->myGrid, tmp_complete, &InnerGrid, &rhs_inner);
	OperationMatrix* L_inner = InnerGrid->createOperationLaplace();
	alpha_inner = new DataVector(rhs_inner->getSize());
	alpha_inner->setAll(0.0);
	rhs_inner->mult(-1.0);
	myCG->solve(*L_inner, *alpha_inner, *rhs_inner, true, verbose, 0.0);

	myGridConverter->updateBoundaryCoefs(alpha, *alpha_inner);

	delete myCG;
	delete L_complete;
	delete L_inner;
	delete rhs_inner;
	delete alpha_inner;
	delete myGridConverter;
	delete myBoundaryUpdate;
	delete InnerGrid;
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
