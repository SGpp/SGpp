/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#define GUNPLOT_RESOLUTION 51
#define SOLUTION_FRAMES 100

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>

#include "sgpp.hpp"

/**
 * Calls the writeHelp method in the BlackScholesSolver Object
 * after creating a screen.
 */
void writeHelp()
{
	sg::LaserHeatEquationSolver* myHESolver = new sg::LaserHeatEquationSolver(1.0, 0.4, 5);

	myHESolver->initScreen();

	delete myHESolver;

	std::stringstream mySStream;

//	mySStream << "Some instructions for the use of Laser Heat Equation Solver:" << std::endl;
//	mySStream << "------------------------------------------------------------" << std::endl << std::endl;
//	mySStream << "Available execution modes are:" << std::endl;



	std::cout << mySStream.str() << std::endl;
}

void testLaserHeatEquation(size_t dim, size_t level, double bound_left, double bound_right, double a,
						double T, double dt, std::string ODESolver,
						double cg_eps, size_t cg_its)
{
	size_t timesteps = (size_t)(T/dt);

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = bound_left;
		myBoundaries[i].rightBoundary = bound_right;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::LaserHeatEquationSolver* myHESolver = new sg::LaserHeatEquationSolver(1.0, 0.04, 5);
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// init Screen Object
	myHESolver->initScreen();

	// Construct a grid
	myHESolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector (start solution)
	DataVector* alpha = new DataVector(myHESolver->getNumberGridPoints());
	myHESolver->refineInitialGridWithLaserHeat(*alpha, 4.0, 5);

	// Print the initial heat function into a gnuplot file
	if (dim < 3)
	{
		myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "laserheatStart.gnuplot");
	}

	// set heat coefficient
	myHESolver->setHeatCoefficient(a);

	// Start solving the Heat Equation
	myHESolver->solveImplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, true, std::max(timesteps/SOLUTION_FRAMES,(size_t)1));

	// Print the solved Heat Equation into a gnuplot file
	if (dim < 3)
	{
		myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "laserheatSolved.gnuplot");
	}

	delete myHESolver;
	delete myBoundingBox;
	delete alpha;
}

int main(int argc, char *argv[])
{
	std::string option;

//	if (argc == 1)
//	{
//		writeHelp();
//		return 0;
//	}

//	option.assign(argv[1]);

//	if (option == "HeatEquation")
//	{
//		if (argc != 13)
//		{
//			writeHelp();
//			return 0;
//		}

		size_t dim;
		size_t level;
		double bound_left;
		double bound_right;
		double a;
		std::string initFunc;
		double T;
		double dt;
		std::string ODESolver;
		double cg_eps;
		size_t cg_its;

		dim = 2;
		level = 6;
		bound_left = 0.0;
		bound_right = 1.0;
		a = 10.0;
		T = 1.0;
		dt = 0.001;
		ODESolver = "ImEul";
		cg_eps = 0.00001;
		cg_its = 400;

		testLaserHeatEquation(dim, level, bound_left, bound_right, a, T, dt, ODESolver, cg_eps, cg_its);
//	}
}
