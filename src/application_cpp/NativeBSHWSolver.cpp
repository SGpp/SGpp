/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi (qic@in.tum.de)

#include "sgpp.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iomanip>

#define CRNIC_IMEUL_STEPS 3

int readStochasticData(std::string tFile, size_t numAssests, DataVector& mu, DataVector& sigmabs, DataMatrix& rho)
{
	std::fstream file;
	double cur_mu;
	double cur_sigmabs;
	double cur_rho;

	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot read file: " << tFile << std::endl;
		return -1;
	}

	for (size_t i = 0; i < numAssests; i++)
	{
		file >> cur_mu;
		file >> cur_sigmabs;
		mu.set(i, cur_mu);
		sigmabs.set(i, cur_sigmabs);
		for (size_t j = 0; j < numAssests; j++)
		{
			file >> cur_rho;
			rho.set(i,j, cur_rho);
		}
	}

	file.close();

	return 0;
}

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
 * Combine Hull White solver and Black Scholes solver test with n assets (ND Sparse Grid) European call option
 *
 * @param l the number of levels used in the Sparse Grid
 * @param the stochastic data (theta, sigmahw, sigmabs, a)
 * @param the grid's bounding box - domain boundary(min,max)
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 */
void testBSHW(size_t d,size_t l, double theta, double sigmahw, double a, std::string fileStoch, std::string fileBound, std::string payoffType,
		size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, double T,double dStrike, bool isLogSolve)
{
	    size_t dim = d;
		size_t level = l;
		size_t timesteps = timeSt;
		double stepsize = dt;
		size_t CGiterations = CGIt;
		double CGepsilon = CGeps;
		DataVector mu(1);
		DataVector sigmabs(1);
		DataMatrix rho(1,1);

		if (readStochasticData(fileStoch, 1, mu, sigmabs, rho) != 0)
			{
				return;
			}

		sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
		if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
		{
			return;
		}

		sg::BlackScholesSolver* myBSSolver;
			if (isLogSolve == true)
			{
				myBSSolver = new sg::BlackScholesSolver(true);
			}
			else
			{
				myBSSolver = new sg::BlackScholesSolver(false);
			}

		sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
		delete[] myBoundaries;

		// init Screen Object
		myBSSolver->initScreen();

		// Construct a grid
		myBSSolver->constructGrid(*myBoundingBox, level);

		// init the basis functions' coefficient vector
		DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

		std::cout << "Grid has " << level << " Levels" << std::endl;
		std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
		std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

		// Init the grid with on payoff function
		myBSSolver->initGridWithPayoffBSHW(*alpha, dStrike, payoffType);

		// Gridpoints @Money
	//	std::cout << "Gridpoints @Money: " << myBSSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;

		// Print the payoff function into a gnuplot file


			myBSSolver->printGrid(*alpha, 20, "payoffBSHW.gnuplot");
			myBSSolver->printSparseGrid(*alpha, "payoffBSHW_surplus.grid.gnuplot", true);
			myBSSolver->printSparseGrid(*alpha, "payoffBSHW_nodal.grid.gnuplot", false);

		//	sg::HullWhiteSolver* myHWSolver = new sg::HullWhiteSolver();

			std::vector<size_t> algoDims(1);
			myBSSolver->setAlgorithmicDimensions(algoDims);

		// Set stochastic data
			myBSSolver->setStochasticData(mu, sigmabs, rho, 0.04);

		//myBSSolver->setStochasticData(0, sigmabs, 0, 0.04);

		// Start solving the Black Scholes Equation
		if (Solver == "ExEul")
		{
			myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
		}
		else if (Solver == "ImEul")
		{
			myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
		}
		else if (Solver == "CrNic")
		{
			myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
		}
		else
		{
			std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
		}

		if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
		{

				// Print the solved Black Scholes Equation into a gnuplot file
				myBSSolver->printGrid(*alpha, 20, "solvedBSHW.gnuplot");
				myBSSolver->printSparseGrid(*alpha, "solvedBSHW_surplus.grid.gnuplot", true);
				myBSSolver->printSparseGrid(*alpha, "solvedBSHW_nodal.grid.gnuplot", false);


		}

		// Test call @ the money
		/*std::vector<double> point;
		for (size_t i = 0; i < d; i++)
		{
			point.push_back(1.0);
		}
		std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;
		*/

		delete alpha;
		delete myBSSolver;
		delete myBoundingBox;

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

	mySStream << "Some instructions for the use of combing Hull White and Black Scholes Solver:" << std::endl;
	mySStream << "------------------------------------------------------" << std::endl << std::endl;
	mySStream << "Available execution modes are:" << std::endl;
	mySStream << "  solveND             Solves an European Call/Put option" << std::endl;
	mySStream << "                      for N assets on a regular sparse grid" << std::endl << std::endl;

	mySStream << "Execution modes descriptions:" << std::endl;
	mySStream << "-----------------------------------------------------" << std::endl;
	mySStream << "solveND" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	value of Theta: theta function" << std::endl;
	mySStream << "	value of sigmahw: sigma value-determine overall level of volatility for hull white" << std::endl;
    mySStream << "	value of a: a" << std::endl;
    mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
    mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	Strike: the strike" << std::endl;
	mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords" << std::endl;

	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "2 5 0.006 0.01 0.1 stoch.data  bound.data " << " std_euro_call "<< "1.0 " << "0.01 "<< "400 " << "0.000001 " << "ExEul " <<"0.3 cart"<<std::endl;
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
		if (argc != 17)
		{
			writeHelp();
		}
		else
		{
			std::string solver;
			std::string payoff;
			double theta;
			double sigmahw;
			double a;
			double dStrike;
			std::string fileStoch;
			std::string fileBound;


			theta = atof(argv[4]);
			sigmahw = atof(argv[5]);
			a = atof(argv[6]);
			fileStoch.assign(argv[7]);
			fileBound.assign(argv[8]);
			payoff.assign(argv[9]);
			solver.assign(argv[14]);
			dStrike = atof(argv[15]);
			std::string coordsType;
			bool coords = false;
			coordsType.assign(argv[16]);
		    if (coordsType == "cart")
			    {
					coords = false;
				}
					else if (coordsType == "log")
				{
					coords = true;
				}
				    else
				{
					std::cout << "Unsupported coordinate option! cart or log are supported!" << std::endl;
					std::cout << std::endl << std::endl;
					writeHelp();
				}
			//testHullWhite(size_t l, double theta, double signma, double a, std::string fileBound, std::string payoffType,
				//	size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, double t, double T)
			testBSHW(atoi(argv[2]),atoi(argv[3]), theta, sigmahw, a, fileStoch, fileBound, payoff, (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[12]), atof(argv[13]), solver,atof(argv[10]),dStrike,coords);
		}
	}

	else
	{
		writeHelp();
	}

	return 0;
}
