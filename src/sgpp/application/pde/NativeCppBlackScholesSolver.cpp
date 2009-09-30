/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "sgpp.hpp"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>

/**
 * reads the values of mu, sigma and rho of all assets from
 * a file and stores them into three separated DataVectors
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssests the of Assets stored in the file
 * @param mu DataVector for the exspected values
 * @param sigma DataVector for standard deviation
 * @param rho DataVector for the correlations
 */
void readStochasticData(std::string tFile, size_t numAssests, DataVector& mu, DataVector& sigma, DataVector& rho)
{
	std::fstream file;
	double cur_mu;
	double cur_sigma;
	double cur_rho;

	file.open(tFile.c_str());

	for (size_t i = 0; i < numAssests; i++)
	{
		file >> cur_mu;
		file >> cur_sigma;
		mu.set(i, cur_mu);
		sigma.set(i, cur_sigma);
		for (size_t j = 0; j < numAssests; j++)
		{
			file >> cur_rho;
			rho.set((i*numAssests)+j, cur_rho);
		}
	}

	file.close();
}

/**
 * reads the values of the Bounding Box
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssests the of Assets stored in the file
 * @param BoundaryArray Pointer to the Bounding Box array
 */
void readBoudingBoxData(std::string tFile, size_t numAssests, sg::DimensionBoundary* BoundaryArray)
{
	std::fstream file;
	double cur_right;
	double cur_left;

	file.open(tFile.c_str());

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
}

/**
 * Do a Black Scholes solver test with one asset (1D Sparse Grid) call option
 *
 * @param l the number of levels used in the Sparse Grid
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param strike1 the strike of the call option
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param animation set this to true if you want to create several pictures during solving in order to create an animation
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 */
void testOneUnderlying(size_t l, std::string fileStoch, std::string fileBound, double strike1, double riskfree, size_t timeSt,
						double dt, size_t CGIt, double CGeps, bool animation, std::string Solver)
{
	size_t dim = 1;
	size_t level = l;
	double* strike = new double[dim];
	strike[0] = strike1;

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataVector rho(dim);

	double r = riskfree;

	readStochasticData(fileStoch, dim, mu, sigma, rho);

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	readBoudingBoxData(fileBound, dim, myBoundaries);

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);

	// init Screen Object
	myBSSolver->initScreen();

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	// Init the grid with on payoff function
	myBSSolver->initGridWithPayoff(*alpha, strike);

	// Print the payoff function into a gnuplot file
	myBSSolver->printGrid(*alpha, 50, "payoff.gnuplot");

	// Set stochastic data
	myBSSolver->setStochasticData(mu, sigma, rho, r);

	// Start solving the Black Scholes Equation
	if (Solver == "ExEul")
	{
		myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, animation, 20);
	}
	else if (Solver == "ImEul")
	{
		myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, animation, 20);
	}
	else if (Solver == "CrNic")
	{
		myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);
	}
	else
	{
		std::cout << "!!!! Use have chosen an unsupoorted solver type !!!!" << std::endl;
	}

	if (Solver == "ImEul" || Solver == "ExEul" || Solver == "CrNic")
	{
		// Print the solved Black Scholes Equation into a gnuplot file
		myBSSolver->printGrid(*alpha, 50, "solvedBS.gnuplot");

		// Do analytic test
		std::vector< std::pair<double, double> >premium;
		double t = (((double)timesteps)*stepsize);
		myBSSolver->solve1DAnalytic(premium, myBoundaries[0].rightBoundary, myBoundaries[0].rightBoundary/50.0, strike[0], t);
		myBSSolver->print1DAnalytic(premium, "analyticBS.gnuplot");
	}

	delete[] myBoundaries;
	delete myBSSolver;
	delete myBoundingBox;
	delete alpha;
}

/**
 * Do a Black Scholes solver test with one asset (1D Sparse Grid) call option
 *
 * @param l the number of levels used in the Sparse Grid
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param strike1 the strike of the call option, first asset
 * @param strike2 the strike of the call option, second asset
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param animation set this to true if you want to create several pictures during solving in order to create an animation
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 */
void testTwoUnderlyings(size_t l, std::string fileStoch, std::string fileBound, double strike1, double strike2, double riskfree, size_t timeSt,
		double dt, size_t CGIt, double CGeps, bool animation, std::string Solver)
{
	size_t dim = 2;
	size_t level = l;
	double* strike = new double[dim];
	strike[0] = strike1;
	strike[1] = strike2;

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataVector rho(dim,dim);

	double r = riskfree;

	readStochasticData(fileStoch, dim, mu, sigma, rho);

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	readBoudingBoxData(fileBound, dim, myBoundaries);

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// init Screen Object
	myBSSolver->initScreen();

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	// Init the grid with on payoff function
	myBSSolver->initGridWithPayoff(*alpha, strike);

	// Print the payoff function into a gnuplot file
	myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");

	// Set stochastic data
	myBSSolver->setStochasticData(mu, sigma, rho, r);

	// Start solving the Black Scholes Equation
	if (Solver == "ExEul")
	{
		myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, animation, 20);
	}
	else if (Solver == "ImEul")
	{
		myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, animation, 20);
	}
	else if (Solver == "CrNic")
	{
		myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);
	}
	else
	{
		std::cout << "!!!! Use have chosen an unsupoorted solver type !!!!" << std::endl;
	}

	if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
	{
		// Print the solved Black Scholes Equation into a gnuplot file
		myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
	}

	delete myBSSolver;
	delete myBoundingBox;
	delete alpha;
}

/**
 * solves a predefined, in the format of University Bonn, grid.
 *
 * @param fileIn the file the contains the grid that should be solved
 * @param fileOut the file the contains the solution grid, written when finished
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 */
void solveBonn(std::string fileIn, std::string fileOut, std::string fileStoch, double riskfree, size_t timeSt,
		double dt, size_t CGIt, double CGeps, std::string Solver)
{
	size_t dim;
	bool hier;

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	double r = riskfree;

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	DataVector* alpha = new DataVector(0);

	// init Screen Object
	myBSSolver->initScreen();

	// Construct a grid, read it from Bonn's format
	myBSSolver->constructGrid(fileIn, *alpha, hier);
	dim = myBSSolver->getNumberDimensions();

	// read stochastic data
	DataVector mu(dim);
	DataVector sigma(dim);
	DataVector rho(dim, dim);
	readStochasticData(fileStoch, dim, mu, sigma, rho);

	// Print the payoff function into a gnuplot file
	//myBSSolver->printGrid(*alpha, 50, "payoff.gnuplot");

	// Set stochastic data
	myBSSolver->setStochasticData(mu, sigma, rho, r);

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
		myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha);
	}
	else
	{
		std::cout << "!!!! Use have chosen an unsupoorted solver type !!!!" << std::endl;
	}

	if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
	{
		// Print the solved Black Scholes Equation into a gnuplot file
		//myBSSolver->printGrid(*alpha, 50, "solvedBS.gnuplot");

		// export the grid, store it to Bonn's format
		myBSSolver->storeGrid(fileOut, *alpha, hier);
	}

	delete myBSSolver;
	delete alpha;
}

/**
 * Calls the writeHelp method in the BlackScholesSolver Object
 * after creating a screen.
 */
void writeHelp()
{
	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();

	// init Screen Object
	myBSSolver->initScreen();

	// print Help
	myBSSolver->writeHelp();

	delete myBSSolver;
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

	if (option == "test1D")
	{
		if (argc != 13)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			std::string ani;
			std::string solver;
			bool animation;

			fileStoch.assign(argv[4]);
			fileBound.assign(argv[3]);
			solver.assign(argv[9]);
			ani.assign(argv[12]);
			if (ani == "animation")
			{
				animation = true;
			}
			else
			{
				animation = false;
			}
			testOneUnderlying(atoi(argv[2]), fileStoch, fileBound, atof(argv[5]), atof(argv[6]), (size_t)(atof(argv[7])/atof(argv[8])), atof(argv[8]), atoi(argv[10]), atof(argv[11]), animation, solver);
		}
	}
	else if (option == "test2D")
	{
		if (argc != 14)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			std::string ani;
			std::string solver;
			bool animation;

			fileStoch.assign(argv[4]);
			fileBound.assign(argv[3]);
			solver.assign(argv[10]);
			ani.assign(argv[13]);
			if (ani == "animation")
			{
				animation = true;
			}
			else
			{
				animation = false;
			}
			testTwoUnderlyings(atoi(argv[2]), fileStoch, fileBound, atof(argv[5]), atof(argv[6]), atof(argv[7]), (size_t)(atof(argv[8])/atof(argv[9])), atof(argv[9]), atoi(argv[11]), atof(argv[12]), animation, solver);
		}
	}
	else if (option == "solveBonn")
	{
		if (argc != 11)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileIn;
			std::string fileOut;
			std::string solver;

			fileIn.assign(argv[2]);
			fileOut.assign(argv[3]);
			fileStoch.assign(argv[4]);
			solver.assign(argv[8]);

			solveBonn(fileIn, fileOut, fileStoch, atof(argv[5]), (size_t)(atof(argv[6])/atof(argv[7])), atof(argv[7]), atoi(argv[9]), atof(argv[10]), solver);
		}
	}
	else
	{
		writeHelp();
	}

	return 0;
}
