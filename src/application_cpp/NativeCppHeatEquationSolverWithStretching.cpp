/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

/// default number of Implicit Euler steps when using Crank Nicolson
#define CRNIC_IMEUL_STEPS 3

#include <time.h>
#include "sgpp.hpp"

//using namespace sg;

void testHeatEquation()
{
	size_t dim = 2;
	size_t level = 8;

	size_t timesteps = 100;
	double stepsize = 0.1;
	size_t CGiterations = 8000;
	double CGepsilon = 0.000001;

	double a = 1.0;

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];

	// set the bounding box
	for (size_t i = 0; i < dim; i++)
	{
		myBoundaries[i].leftBoundary = 0.001;
		myBoundaries[i].rightBoundary = 3.0;
		myBoundaries[i].bDirichletLeft = true;
		myBoundaries[i].bDirichletRight = true;
	}

	sg::Stretching1D* stretching1Ds = new sg::Stretching1D[dim];
		string s0("id");
		string s1("log");
		string s2("sinh");

		for(size_t j=0;j<dim;j++){
			stretching1Ds[j].type.assign(s1);
			stretching1Ds[j].x_0 = 0.5;
			stretching1Ds[j].xsi=1.0;
		}

	sg::HeatEquationSolverWithStretching* myHESolver = new sg::HeatEquationSolverWithStretching();
	sg::Stretching* myStretching = new sg::Stretching(dim, myBoundaries, stretching1Ds);
	delete[] myBoundaries;
	delete[] stretching1Ds;
//	myStretching->printLookupTable();
	// init Screen Object
	myHESolver->initScreen();

	// Construct a grid
	myHESolver->constructGrid(*myStretching, level);


	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myHESolver->getNumberGridPoints());

//	myHESolver->initGridWithSingleHeat(*alpha, 100.0);
	myHESolver->initGridWithSmoothHeat(*alpha, 3.0, 1.5, 5.0);
//	myHESolver->initGridWithConstantHeat(*alpha, 4.0);

//	cout<< (*alpha).toString() << endl;


	// Print the initial heat function into a gnuplot file
	myHESolver->printGrid(*alpha, 50, "heatStart.gnuplot");

	// set heat coefficient
	myHESolver->setHeatCoefficient(a);

	// Start solving the Heat Equation
	//myHESolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, true, true, 50);
	myHESolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, true, false, 50);
//	myHESolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);

	// Print the solved Heat Equation into a gnuplot file
	myHESolver->printGrid(*alpha, 50, "solvedHeat.gnuplot");

	////Problem with deconstructing grid, uncomment when that is solved.
//	delete myHESolver;
	delete myStretching;
	delete alpha;
}


int main(int argc, char *argv[])
{
	clock_t start =clock();
	testHeatEquation();
	clock_t end =clock();
	std::cout<<"runtime: "<<(double)(end-start)/CLOCKS_PER_SEC<<"seconds"<<std::endl;
}
