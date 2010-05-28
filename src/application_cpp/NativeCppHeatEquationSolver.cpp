/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

void testHeatEquation()
{
	size_t dim = 2;
	size_t level = 7;

	size_t timesteps = 100;
	double stepsize = 0.1;
	size_t CGiterations = 8000;
	double CGepsilon = 0.000001;

	double a = 1.0;

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = 0.0;
		myBoundaries[i].rightBoundary = 3.0;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::HeatEquationSolver* myHESolver = new sg::HeatEquationSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// init Screen Object
	myHESolver->initScreen();

	// Construct a grid
	myHESolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myHESolver->getNumberGridPoints());

	//myHESolver->initGridWithSingleHeat(*alpha, 100.0);
	myHESolver->initGridWithSmoothHeat(*alpha, 3.0, 1.5, 5.0);
	//myHESolver->initGridWithConstantHeat(*alpha, 4.0);

	// Print the initial heat function into a gnuplot file
	myHESolver->printGrid(*alpha, 50, "heatStart.gnuplot");

	// set heat coefficient
	myHESolver->setHeatCoefficient(a);

	// Start solving the Heat Equation
	//myHESolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, true, true, 50);
	myHESolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, true, false, 50);
	//myHESolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);

	// Print the solved Heat Equation into a gnuplot file
	myHESolver->printGrid(*alpha, 50, "solvedHeat.gnuplot");

	delete myHESolver;
	delete myBoundingBox;
	delete alpha;
}


int main(int argc, char *argv[])
{
	testHeatEquation();
}
