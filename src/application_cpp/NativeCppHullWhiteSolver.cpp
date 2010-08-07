/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi (qic@in.tum.de)

#include "sgpp.hpp"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iomanip>

#define CRNIC_IMEUL_STEPS 3

int readBoudingBoxData(std::string tFile, size_t numAssests, sg::DimensionBoundary* BoundaryArray)
{
	std::fstream file;
	double cur_right;
	double cur_left;

	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot read file: " << tFile << std::endl;
		return -1;
	}

	for (size_t i = 0; i < numAssests; i++)
	{
		file >> cur_left;
		file >> cur_right;

		BoundaryArray[i].leftBoundary = cur_left;
		BoundaryArray[i].rightBoundary = cur_right;
		BoundaryArray[i].bDirichletLeft = true;
		BoundaryArray[i].bDirichletRight = true;
	}

	file.close();

	return 0;
}

/**
 * Do a Hull White solver test with n assets (ND Sparse Grid) European call option
 *
 * @param l the number of levels used in the Sparse Grid
 * @param the stochastic data (theta, sigma, a)
 * @param the grid's bounding box - domain boundary(min,max)
 * @param r the interest rate
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 */
void testHullWhite(size_t l, double theta, double sigma, double a, std::string fileBound, std::string payoffType,
		size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, double t, double T,double dStrike)
{
		size_t level = l;
		size_t timesteps = timeSt;
		double stepsize = dt;
		size_t CGiterations = CGIt;
		double CGepsilon = CGeps;


		sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[1];
		if (readBoudingBoxData(fileBound, 1, myBoundaries) != 0)
		{
			return;
		}

		sg::HullWhiteSolver* myHWSolver = new sg::HullWhiteSolver();
		sg::BoundingBox* myBoundingBox = new sg::BoundingBox(1, myBoundaries);
		delete[] myBoundaries;

		// init Screen Object
		myHWSolver->initScreen();

		// Construct a grid
		myHWSolver->constructGrid(*myBoundingBox, level);

		// init the basis functions' coefficient vector
		DataVector* alpha = new DataVector(myHWSolver->getNumberGridPoints());

		std::cout << "Grid has " << level << " Levels" << std::endl;
		std::cout << "Initial Grid size: " << myHWSolver->getNumberGridPoints() << std::endl;
		std::cout << "Initial Grid size (inner): " << myHWSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

		// Init the grid with on payoff function
		myHWSolver->initGridWithPayoff(*alpha, dStrike, payoffType, sigma, a, t, T);

		// Gridpoints @Money
	//	std::cout << "Gridpoints @Money: " << myBSSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;

		// Print the payoff function into a gnuplot file


			myHWSolver->printGrid(*alpha, 20, "payoffHW.gnuplot");
			myHWSolver->printSparseGrid(*alpha, "payoffHW_surplus.grid.gnuplot", true);
			myHWSolver->printSparseGrid(*alpha, "payoffHW_nodal.grid.gnuplot", false);


		// Set stochastic data
		//myBSSolver->setStochasticData(mu, sigma, rho, r);

		// Start solving the Black Scholes Equation
		if (Solver == "ExEul")
		{
			myHWSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
		}
		else if (Solver == "ImEul")
		{
			myHWSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
		}
		else if (Solver == "CrNic")
		{
			myHWSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
		}
		/*else if (Solver == "AdBas")
		{
			myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
		}
		else if (Solver == "VaTim")
		{
			myBSSolver->solveVarTimestep(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
		}*/
		else
		{
			std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
		}

		if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
		{

				// Print the solved Black Scholes Equation into a gnuplot file
				myHWSolver->printGrid(*alpha, 20, "solvedHW.gnuplot");
				myHWSolver->printSparseGrid(*alpha, "solvedHW_surplus.grid.gnuplot", true);
				myHWSolver->printSparseGrid(*alpha, "solvedHW_nodal.grid.gnuplot", false);


		}

		// Test call @ the money
		/*std::vector<double> point;
		for (size_t i = 0; i < d; i++)
		{
			point.push_back(1.0);
		}
		std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;

		delete alpha;
		delete myBSSolver;
		delete myBoundingBox;
*/
	//std::cout << "Nothing has been implemented so far " << std::endl;


}

/**
 * Calls the writeHelp method in the HullWhiteSolver Object
 * after creating a screen.
 */

void writeHelp()
{
	/*sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();

	myBSSolver->initScreen();

	delete myBSSolver;*/

	std::stringstream mySStream;

	mySStream << "Some instructions for the use of Hull White Solver:" << std::endl;
	mySStream << "------------------------------------------------------" << std::endl << std::endl;
	mySStream << "Available execution modes are:" << std::endl;
	mySStream << "  solveND             Solves an European Call/Put option" << std::endl;
	mySStream << "                      for N assets on a regular sparse grid" << std::endl << std::endl;

	mySStream << "Execution modes descriptions:" << std::endl;
	mySStream << "-----------------------------------------------------" << std::endl;
	mySStream << "solveND" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	value of Theta: theta function" << std::endl;
	mySStream << "	value of sigma: sigma value-determine overall level of volatility" << std::endl;
    mySStream << "	value of a: a" << std::endl;
    mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	//mySStream << "	r: the interest rate" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	t: the current time" << std::endl;
	mySStream << "	Strike: the strike" << std::endl;

	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "5 0 0.01 0.1 bound.data " << " std_euro_call "<< "1.0 " << "0.01 "<< "400 " << "0.000001 " << "ImEul " << "0.5 " <<"1.0 "<<std::endl;
	mySStream << std::endl;

	mySStream << std::endl << std::endl;
	std::cout << mySStream.str() << std::endl;
}

/**
 * main routine of the application, do some first cli
 * correction test and branches to right solver configuration
 *
 * @param argc contains the number of cli arguments
 * @param argv contains the cli arguments as C-Strings
 */
int main(int argc, char *argv[])
{
	std::string option;

	if (argc == 1)
	{
		writeHelp();

		return 0;
	}

	option.assign(argv[1]);

	if (option == "solveND")
	{
		if (argc != 15)
		{
			writeHelp();
		}
		else
		{
			std::string solver;
			std::string payoff;
			double theta;
			double sigma;
			double a;
			double dStrike;
			std::string fileBound;


			theta = atof(argv[3]);
			sigma = atof(argv[4]);
			a = atof(argv[5]);
			fileBound.assign(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[12]);
			dStrike = atof(argv[14]);
			//testHullWhite(size_t l, double theta, double signma, double a, std::string fileBound, std::string payoffType,
				//	size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, double t, double T)
			testHullWhite(atoi(argv[2]), theta, sigma, a, fileBound, payoff, (size_t)(atof(argv[8])/atof(argv[9])), atof(argv[9]), atoi(argv[10]), atof(argv[11]), solver,atof(argv[13]),atof(argv[8]),dStrike);
		}
	}

	else
	{
		writeHelp();
	}

	return 0;
}
