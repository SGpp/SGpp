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
void testHullWhite(size_t l, double theta, double signma, double a, double min, double max, double r, std::string payoffType,
		size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver)
{

	std::cout << "Nothing has been implemented so far " << std::endl;


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
	mySStream << "	value of min: the left range of the boundary" << std::endl;
	mySStream << "	value of max: the right range of the boundary" << std::endl;
	mySStream << "	r: the interest rate" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "5 0 0.01 0.1 0 2.5 0.04 " << " std_euro_call "<< "1.0 " << "0.01 "<< "400 " << "0.000001 " << "ImEul" << std::endl;
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
			double min;
			double max;
			double r;

			theta = atof(argv[3]);
			sigma = atof(argv[4]);
			a = atof(argv[5]);
			min = atof(argv[6]);
			max = atof(argv[7]);
			r = atof(argv[8]);
			payoff.assign(argv[9]);
			solver.assign(argv[14]);

			testHullWhite(atoi(argv[2]), theta, sigma, a, min, max, r, payoff, (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[12]), atof(argv[13]), solver);
		}
	}

	else
	{
		writeHelp();
	}

	return 0;
}
