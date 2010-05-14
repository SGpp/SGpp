/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
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
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readStochasticData(std::string tFile, size_t numAssests, DataVector& mu, DataVector& sigma, DataVector& rho)
{
	std::fstream file;
	double cur_mu;
	double cur_sigma;
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

	return 0;
}

/**
 * reads the values of the Bounding Box
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssests the of Assets stored in the file
 * @param BoundaryArray Pointer to the Bounding Box array
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
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
 * reads the strikes of all assets from
 * a file and stores them into an array
 *
 * @param tFile the file that contains the strikes
 * @param numAssests the of Assets stored in the file
 * @param strikes array into which the stikes should be stored
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readStrikes(std::string tFile, size_t numAssests, double* strikes)
{
	std::fstream file;
	double cur_strike;

	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot read file: " << tFile << std::endl;
		return -1;
	}

	for (size_t i = 0; i < numAssests; i++)
	{
		file >> cur_strike;
		strikes[i] = cur_strike;
	}

	file.close();

	return 0;
}

/**
 * reads the analyze configuration from a file
 *
 * @param tFile the file that contains the analyze data
 * @param numAssests the number of assets
 * @param percent variable to store size of cuboid in every dimension
 * @param points variable to store the number of points in every dimension
 * @param center vector to store the center of the evaluation cuboid
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readAnalyzeData(std::string tFile, size_t numAssests, double& percent, size_t& points, std::vector<double>& center)
{
	std::fstream file;
	double cur_coord;

	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot read file: " << tFile << std::endl;
		return -1;
	}

	center.empty();
	for (size_t i = 0; i < numAssests; i++)
	{
		file >> cur_coord;
		center.push_back(cur_coord);
	}

	file >> percent;
	file >> points;

	file.close();

	return 0;
}


/**
 * Do a Black Scholes solver test with one asset (1D Sparse Grid) European call option
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

	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
	{
		return;
	}

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
	{
		return;
	}

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);

	// init Screen Object
	myBSSolver->initScreen();

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	// Init the grid with on payoff function
	myBSSolver->initGridWithEuroCallPayoff(*alpha, strike, "avg");

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
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
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
 * Do a Black Scholes solver test with one asset (2D Sparse Grid) European call option
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

	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
	{
		return;
	}

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
	{
		return;
	}

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
	myBSSolver->initGridWithEuroCallPayoff(*alpha, strike, "max");

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
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
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
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call option
 *
 * @param d the number of dimensions used in the Sparse Grid
 * @param l the number of levels used in the Sparse Grid
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param fileStrike filename of the file that contains the assets' strikes
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 */
void testNUnderlyings(size_t d, size_t l, std::string fileStoch, std::string fileBound, std::string fileStrike, std::string payoffType,
		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver)
{
	size_t dim = d;
	size_t level = l;
	double* strike = new double[dim];

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataVector rho(dim,dim);

	double r = riskfree;

	if (readStrikes(fileStrike, dim, strike) != 0)
	{
		return;
	}

	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
	{
		return;
	}

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
	{
		return;
	}

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
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
	myBSSolver->initGridWithEuroCallPayoff(*alpha, strike, payoffType);

	// Print the payoff function into a gnuplot file
	if (dim < 3)
	{
		myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");
		myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
		myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
	}
	else
	{
		myBSSolver->storeGridBonn("payoff_Nd.bonn", *alpha, true);
	}

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
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
	{
		if (dim < 3)
		{
			// Print the solved Black Scholes Equation into a gnuplot file
			myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
			myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
			myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
		}
		else
		{
			myBSSolver->storeGridBonn("solvedBS_Nd.bonn", *alpha, true);
		}
	}

	// Test call @ the money
	if (payoffType == "avgM")
	{
		std::vector<double> point;
		for (size_t i = 0; i < d; i++)
		{
			point.push_back(1.0);
		}
		std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;
	}

	delete alpha;
	delete myBSSolver;
	delete myBoundingBox;
}

/**
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call option
 *
 * @param d the number of dimensions used in the Sparse Grid
 * @param start_l the number of levels used in the Sparse Grid (first test)
 * @param end_l the number of level used in the Sparse Grid (last test)
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param fileStrike filename of the file that contains the assets' strikes
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 */
void testNUnderlyingsAnalyze(size_t d, size_t start_l, size_t end_l, std::string fileStoch, std::string fileBound, std::string fileStrike, std::string payoffType,
		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, std::string fileAnalyze)
{
	size_t dim = d;
	double* strike = new double[dim];

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataVector rho(dim,dim);

	DataVector EvalPoints(1, d);

	double r = riskfree;

	std::vector<DataVector> results;

	if (readStrikes(fileStrike, dim, strike) != 0)
	{
		return;
	}

	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
	{
		return;
	}

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
	{
		return;
	}

	double cuboidSize = 0.0;
	size_t points = 0;
 	std::vector<double> center;
	if (readAnalyzeData(fileAnalyze, dim, cuboidSize, points, center) != 0)
	{
		return;
	}

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// init Screen Object
	myBSSolver->initScreen();

	for (size_t i = start_l; i <= end_l; i++)
	{
		size_t level = i;

		// Construct a grid
		myBSSolver->constructGrid(*myBoundingBox, level);

		// in the first iteration -> calculate the evaluation points
		if (i == start_l)
		{
			myBSSolver->getEvaluationCuboid(EvalPoints, center, cuboidSize, points);

			//std::cout << EvalPoints.toString() << std::endl;
		}

		// init the basis functions' coefficient vector
		DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

		std::cout << "Grid has " << level << " Levels" << std::endl;
		std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
		std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

		if (i > start_l)
			std::cout << std::endl << std::endl << std::endl;

		// Init the grid with on payoff function
		myBSSolver->initGridWithEuroCallPayoff(*alpha, strike, payoffType);

		// Print the payoff function into a gnuplot file
		if (dim < 3)
		{
			myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");
			myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
			myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
		}
		else
		{
			myBSSolver->storeGridBonn("payoff_Nd.bonn", *alpha, true);
		}

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
			std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
		}

		if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
		{
			if (dim < 3)
			{
				// Print the solved Black Scholes Equation into a gnuplot file
				myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
				myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
				myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
			}
			else
			{
				myBSSolver->storeGridBonn("solvedBS_Nd.bonn", *alpha, true);
			}
		}

		// Test call @ the money
		if (payoffType == "avgM")
		{
			std::vector<double> point;
			for (size_t i = 0; i < d; i++)
			{
				point.push_back(1.0);
			}
			std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;
		}

		// Evaluate Cuboid
		DataVector Prices(EvalPoints.getSize());
		myBSSolver->evaluateCuboid(*alpha, Prices, EvalPoints);
		results.push_back(Prices);

		//std::cout << Prices.toString() << std::endl;

		myBSSolver->deleteGrid();
		delete alpha;

		std::cout << std::endl;
	}

	delete myBSSolver;
	delete myBoundingBox;

	std::cout << "=====================================================================" << std::endl;
	std::cout << "=====================================================================" << std::endl << std::endl;
	std::cout << "Calculating norms of relative errors to a grid" << std::endl;
	std::cout << "with " << end_l << " levels and testing-coboid" << std::endl;
	std::cout << "with the center:" << std::endl;
	for (size_t i = 0; i < d; i++)
	{
		std::cout << center[i] << " ";
	}
	std::cout << std::endl << "and " << points << " test-points in a range of " << std::endl;
	std::cout << cuboidSize*200.0 << "% per dimension:" << std::endl << std::endl;

	// Calculate relative errors and some norms
	for (size_t i = 0; i < end_l-start_l; i++)
	{
		DataVector maxLevel(results[end_l-start_l]);
		DataVector relError(results[i]);
		double maxNorm = 0.0;
		double twoNorm = 0.0;

		// calculate relative error
		relError.sub(maxLevel);
		relError.componentwise_div(maxLevel);

		// calculate max. norm of relative error
		maxNorm = relError.max();

		// calculate two norm of relative error
		relError.componentwise_mult(relError);
		twoNorm = relError.sum();
		twoNorm = sqrt(twoNorm);

		// Printing norms
		std::cout << "Level " << i + start_l << ": max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << twoNorm << std::endl;
	}
	std::cout << std::endl << std::endl;
}

/**
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call option, with Intial
 * Grid Refinement
 *
 * @param d the number of dimensions used in the Sparse Grid
 * @param l the number of levels used in the Sparse Grid
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param fileStrike filename of the file that contains the assets' strikes
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 * @param nIterAdaptSteps number of the iterative Grid Refinement that should be executed
 * @param dInitialAdpatDist initial distance from @the money. Is devided in every iteration by the number of the iteration
 */
void testNUnderlyingsAdapt(size_t d, size_t l, std::string fileStoch, std::string fileBound, std::string fileStrike, std::string payoffType,
		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, size_t nIterAdaptSteps, double dInitialAdpatDist)
{
	size_t dim = d;
	size_t level = l;
	double* strike = new double[dim];

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataVector rho(dim,dim);

	double r = riskfree;

	if (readStrikes(fileStrike, dim, strike) != 0)
	{
		return;
	}

	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
	{
		return;
	}

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
	{
		return;
	}

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// init Screen Object
	myBSSolver->initScreen();

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
	std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

	// refine the grid to approximate the singularity in the start solution better
	for (size_t i = 0 ; i < nIterAdaptSteps; i++)
	{
		std::cout << "Refining Grid..." << std::endl;
		myBSSolver->refineInitialGrid(*alpha, strike, payoffType, (dInitialAdpatDist/(static_cast<double>(i+1))));
		std::cout << "Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
	}
	std::cout << std::endl << std::endl;

	// Print the payoff function into a gnuplot file
	if (dim < 3)
	{
		myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");
		myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
		myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
	}
	else
	{
		myBSSolver->storeGridBonn("payoff_Nd.bonn", *alpha, true);
	}

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
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
	{
		if (dim < 3)
		{
			// Print the solved Black Scholes Equation into a gnuplot file
			myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
			myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
			myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
		}
		else
		{
			myBSSolver->storeGridBonn("solvedBS_Nd.bonn", *alpha, true);
		}
	}

	// Test call @ the money
	if (payoffType == "avgM")
	{
		std::vector<double> point;
		for (size_t i = 0; i < d; i++)
		{
			point.push_back(1.0);
		}
		std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;
	}

	delete myBSSolver;
	delete myBoundingBox;
	delete alpha;
}

/**
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call option, with Intial
 * Grid Refinement
 *
 * @param d the number of dimensions used in the Sparse Grid
 * @param l the number of levels used in the Sparse Grid
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param fileStrike filename of the file that contains the assets' strikes
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 * @param refinePercent percantage of points that should be refined before Black Scholes Equation is solved
 */
void testNUnderlyingsAdaptSurplus(size_t d, size_t l, std::string fileStoch, std::string fileBound, std::string fileStrike, std::string payoffType,
		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, double refinePercent)
{
	size_t dim = d;
	size_t level = l;
	double* strike = new double[dim];

	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataVector rho(dim,dim);

	double r = riskfree;

	if (readStrikes(fileStrike, dim, strike) != 0)
	{
		return;
	}

	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
	{
		return;
	}

	sg::DimensionBoundary* myBoundaries = new sg::DimensionBoundary[dim];
	if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
	{
		return;
	}

	sg::BlackScholesSolver* myBSSolver = new sg::BlackScholesSolver();
	sg::BoundingBox* myBoundingBox = new sg::BoundingBox(dim, myBoundaries);
	delete[] myBoundaries;

	// init Screen Object
	myBSSolver->initScreen();

	// Construct a grid
	myBSSolver->constructGrid(*myBoundingBox, level);

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
	std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

	// Init the grid with on payoff function
	myBSSolver->initGridWithEuroCallPayoff(*alpha, strike, payoffType);

	// refine the grid to approximate the singularity in the start solution better
	myBSSolver->refineInitialGridSurplus(*alpha, refinePercent);
	std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
	std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

	// Print the payoff function into a gnuplot file
	if (dim < 3)
	{
		myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");
		myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
		myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
	}
	else
	{
		myBSSolver->storeGridBonn("payoff_Nd.bonn", *alpha, true);
	}

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
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
	{
		if (dim < 3)
		{
			// Print the solved Black Scholes Equation into a gnuplot file
			myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
			myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
			myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
		}
		else
		{
			myBSSolver->storeGridBonn("solvedBS_Nd.bonn", *alpha, true);
		}
	}

	// Test call @ the money
	if (payoffType == "avgM")
	{
		std::vector<double> point;
		for (size_t i = 0; i < d; i++)
		{
			point.push_back(1.0);
		}
		std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;
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
	myBSSolver->constructGridBonn(fileIn, *alpha, hier);
	dim = myBSSolver->getNumberDimensions();

	// read stochastic data
	DataVector mu(dim);
	DataVector sigma(dim);
	DataVector rho(dim, dim);
	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
	{
		return;
	}

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
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic")
	{
		// Print the solved Black Scholes Equation into a gnuplot file
		//myBSSolver->printGrid(*alpha, 50, "solvedBS.gnuplot");

		// export the grid, store it to Bonn's format
		myBSSolver->storeGridBonn(fileOut, *alpha, hier);
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

	myBSSolver->initScreen();

	delete myBSSolver;

	std::stringstream mySStream;

	mySStream << "Some instructions for the use of Black Scholes Solver:" << std::endl;
	mySStream << "------------------------------------------------------" << std::endl << std::endl;
	mySStream << "Available execution modes are:" << std::endl;
	mySStream << "	test1D		Solves a simple 1D example" << std::endl;
	mySStream << "	test2D		Solves a 2D example" << std::endl;
	mySStream << "	solveND		Solves a ND example" << std::endl;
	mySStream << "	solveBonn	Solves an option delivered in Bonn's format" << std::endl << std::endl;

	mySStream << "test1D" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	strike: the strike of the call option" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	[no]animation: generate pictures for an animation" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "5 " << "bound.data stoch.data " << "65.0 " << "0.05 " << "1.0 " << "0.1 ImEul " << "400 " << "0.000001 " << "noanimation" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following output files:" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << "	analyticBS.gnuplot: the analytical solution" << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "test2D" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	strike1: the strike 1 of the call option" << std::endl;
	mySStream << "	strike2: the strike 1 of the call option" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	[no]animation: generate pictures for an animation" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "5 " << "bound.data stoch.data " << "65.0 55.0 " << "0.05 " << "1.0 " << "0.1 ImEul " << "400 " << "0.000001 " << "noanimation" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following output files:" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "solveND" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	file_Strikes: file containing strikes of Europ. Call" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: max, avg/avgM" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 5 " << "bound.data stoch.data strike.data max "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 " << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
	mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
	mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "solveNDanalyze" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level_start: number of levels within the Sparse Grid (start)" << std::endl;
	mySStream << "	level_end: number of levels within the Sparse Grid (end)" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	file_Strikes: file containing strikes of Europ. Call" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: max, avg/avgM" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	file_analyze: file containing the analyzing options" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 2 5 " << "bound.data stoch.data strike.data max "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 anal.data" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
	mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
	mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "solveNDadapt" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	file_Strikes: file containing strikes of Europ. Call" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: max, avg/avgM" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	Adapt-Initial-Refinement: Number of Initial" << std::endl;
	mySStream << "			Refinements" << std::endl;
	mySStream << "	Adapt-Initial-Distance: determines the distance" << std::endl;
	mySStream << "			a grid point must have from @money to" << std::endl;
	mySStream << "			by refined" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 5 " << "bound.data stoch.data strike.data max "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 5 0.5" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
	mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
	mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "solveNDadaptSurplus" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	file_Strikes: file containing strikes of Europ. Call" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: max, avg/avgM" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	Adapt-Refinement Percent: Percent of grid points that should be refined" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 5 " << "bound.data stoch.data strike.data max "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 0.5" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
	mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
	mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;


	mySStream << "solveBonn" << std::endl << "---------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	file_grid_in: file the specifies the unsolved grid" << std::endl;
	mySStream << "	file_grid_out: file that contains the solved grid when finished" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "grid.in grid.out " << "stoch.data " << "0.05 " << "1.0 " << "0.1 ImEul " << "400 " << "0.000001 " << std::endl;

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
	else if (option == "solveND")
	{
		if (argc != 14)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			std::string fileStrike;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[5]);
			fileBound.assign(argv[4]);
			fileStrike.assign(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[11]);

			testNUnderlyings(atoi(argv[2]), atoi(argv[3]), fileStoch, fileBound, fileStrike, payoff, atof(argv[8]), (size_t)(atof(argv[9])/atof(argv[10])), atof(argv[10]), atoi(argv[12]), atof(argv[13]), solver);
		}
	}
	else if (option == "solveNDanalyze")
	{
		if (argc != 16)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			std::string fileStrike;
			std::string fileAnalyze;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[6]);
			fileBound.assign(argv[5]);
			fileStrike.assign(argv[7]);
			fileAnalyze.assign(argv[15]);
			payoff.assign(argv[8]);
			solver.assign(argv[12]);

			testNUnderlyingsAnalyze(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), fileStoch, fileBound, fileStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, fileAnalyze);
		}
	}
	else if (option == "solveNDadapt")
	{
		if (argc != 16)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			std::string fileStrike;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[5]);
			fileBound.assign(argv[4]);
			fileStrike.assign(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[11]);

			testNUnderlyingsAdapt(atoi(argv[2]), atoi(argv[3]), fileStoch, fileBound, fileStrike, payoff, atof(argv[8]), (size_t)(atof(argv[9])/atof(argv[10])), atof(argv[10]), atoi(argv[12]), atof(argv[13]), solver, atoi(argv[14]), atof(argv[15]));
		}
	}
	else if (option == "solveNDadaptSurplus")
	{
		if (argc != 15)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			std::string fileStrike;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[5]);
			fileBound.assign(argv[4]);
			fileStrike.assign(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[11]);

			testNUnderlyingsAdaptSurplus(atoi(argv[2]), atoi(argv[3]), fileStoch, fileBound, fileStrike, payoff, atof(argv[8]), (size_t)(atof(argv[9])/atof(argv[10])), atof(argv[10]), atoi(argv[12]), atof(argv[13]), solver, atoi(argv[14]));
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
