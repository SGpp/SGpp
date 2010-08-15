/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <cmath>

// @todo (heinecke) remove global variables
std::string tFileEvalCuboid = "evalCuboid.data";
std::string tFileEvalCuboidValues = "evalCuboidValues.data";

/// default number of Implicit Euler steps before starting with Crank Nicolson approach
#define CRNIC_IMEUL_STEPS 3
/// default value for epsilon in gridpoints @money
#define DFLT_EPS_AT_MONEY 0.0

/**
 * reads the values of mu, sigma and rho of all assets from
 * a file and stores them into three separated DataVectors
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssests the of Assets stored in the file
 * @param mu DataVector for the exspected values
 * @param sigma DataVector for standard deviation
 * @param rho DataMatrix for the correlations
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readStochasticData(std::string tFile, size_t numAssests, DataVector& mu, DataVector& sigma, DataMatrix& rho)
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
			rho.set(i,j, cur_rho);
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
 * reads a cuboid defined by several points from a file. These points are stored in the
 * cuboid DataMatrix
 *
 * @param cuboid DataMatrix into which the evaluations points are stored
 * @param tFile file that contains the cuboid
 * @param dim the dimensions of cuboid
 */
int readEvalutionCuboid(DataMatrix& cuboid, std::string tFile, size_t dim)
{
	std::fstream file;
	double cur_coord;

	file.open(tFile.c_str());

	if(cuboid.getNcols() != dim)
	{
		std::cout << "Cuboid-definition file doesn't match: " << tFile << std::endl;
		return -1;
	}

	if(!file.is_open())
	{
		std::cout << "Error cannot read file: " << tFile << std::endl;
		return -1;
	}

	size_t i = 0;
	while (!file.eof())
	{
		DataVector line(dim);
		for (size_t d = 0; d < dim; d++)
		{
			file >> cur_coord;
			line.set(d, cur_coord);
		}
		cuboid.resize(i+1);
		cuboid.setRow(i, line);
		i++;
	}

	file.close();

	return 0;
}

/**
 * reads function values (here option prices) from a file
 *
 * @param values DataVector into which the values will be stored
 * @param tFile file from which the values are read
 * @param numValues number of values stored in the file
 */
int readOptionsValues(DataVector& values, std::string tFile)
{
	std::fstream file;
	double cur_value;

	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot read file: " << tFile << std::endl;
		return -1;
	}

	size_t i = 0;
	while (!file.eof())
	{
		file >> cur_value;
		values.resize(i+1);
		values.set(i, cur_value);
		i++;
	}

	file.close();

	return 0;
}

/**
 * Writes a DataMatrix into a file
 *
 * @param data the DataMatrix that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
int writeDataMatrix(DataMatrix& data, std::string tFile)
{
	std::ofstream file;
	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot write file: " << tFile << std::endl;
		return -1;
	}

	for (size_t i = 0; i < data.getNrows(); i++)
	{
		for (size_t j = 0; j < data.getNcols(); j++)
		{
			file << std::scientific << std::setprecision( 16 ) << data.get(i,j) << " ";
		}
		file << std::endl;
	}

	file.close();

	return 0;
}


/**
 * Writes a DataVector into a file
 *
 * @param data the DataVector that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
int writeDataVector(DataVector& data, std::string tFile)
{
	std::ofstream file;
	file.open(tFile.c_str());

	if(!file.is_open())
	{
		std::cout << "Error cannot write file: " << tFile << std::endl;
		return -1;
	}

	for (size_t i = 0; i < data.getSize(); i++)
	{

		file << std::scientific << std::setprecision( 16 ) << data.get(i) << " " << std::endl;
	}

	file.close();

	return 0;
}

/**
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call / put option
 *
 * @param d the number of dimensions used in the Sparse Grid
 * @param l the number of levels used in the Sparse Grid
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param dStrike the strike of the option
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 * @param isLogSolve set this to true if the log-transformed Black Scholes Equation should be solved
 */
void testNUnderlyings(size_t d, size_t l, std::string fileStoch, std::string fileBound, double dStrike, std::string payoffType,
		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, bool isLogSolve)
{
	size_t dim = d;
	size_t level = l;
	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataMatrix rho(dim,dim);

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
	myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);

	// Gridpoints @Money
	std::cout << "Gridpoints @Money: " << myBSSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;

	// Print the payoff function into a gnuplot file
	if (dim < 3)
	{
		if (dim == 2)
		{
			sg::DimensionBoundary* myAreaBoundaries = new sg::DimensionBoundary[dim];

			for (size_t i = 0; i < 2; i++)
			{
				myAreaBoundaries[i].leftBoundary = 0.9;
				myAreaBoundaries[i].rightBoundary = 1.1;
			}
			sg::BoundingBox* myGridArea = new sg::BoundingBox(dim, myAreaBoundaries);

			myBSSolver->printGridDomain(*alpha, 50, *myGridArea, "payoff_area.gnuplot");

			delete[] myAreaBoundaries;
			delete myGridArea;
		}
		myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");
		myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
		myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
		myBSSolver->printPayoffInterpolationError2D(*alpha, "payoff_interpolation_error.grid.gnuplot", 10000, dStrike);
		if (isLogSolve == true)
		{
			myBSSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.gnuplot", true);
			myBSSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.gnuplot", false);
		}
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
		myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
	}
	else if (Solver == "AdBas")
	{
		myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else if (Solver == "VaTim")
	{
		myBSSolver->solveVarTimestep(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else
	{
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (dim < 3)
	{
		// Print the solved Black Scholes Equation into a gnuplot file
		myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
		myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
		myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
		if (isLogSolve == true)
		{
			myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_surplus_cart.grid.gnuplot", true);
			myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_nodal_cart.grid.gnuplot", false);
		}
	}

	// Test call @ the money
	std::vector<double> point;
	for (size_t i = 0; i < d; i++)
	{
		if (isLogSolve == true)
		{
			point.push_back(log(1.0));
		}
		else
		{
			point.push_back(1.0);
		}
	}
	std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;

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
 * @param dStrike the strike of the option
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 * @param isLogSolve set this to true if the log-transformed Black Scholes Equation should be solved
 */
void testNUnderlyingsAnalyze(size_t d, size_t start_l, size_t end_l, std::string fileStoch, std::string fileBound, double dStrike, std::string payoffType,
		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, std::string fileAnalyze, bool isLogSolve)
{
	size_t dim = d;
	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataMatrix rho(dim,dim);

	DataMatrix EvalPoints(1, d);

	double r = riskfree;

	std::vector<DataVector> results;


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
	sg::EvalCuboidGenerator* myEvalCuboidGen = new sg::EvalCuboidGenerator(*myBoundingBox, dim);
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
			myEvalCuboidGen->getEvaluationCuboid(EvalPoints, center, cuboidSize, points);

			writeDataMatrix(EvalPoints, tFileEvalCuboid);

			// If the log-transformed Black Scholes Eqaution is used -> transform Eval-cuboid
			if (isLogSolve == true)
			{
				for (size_t v = 0; v < EvalPoints.getNrows(); v++)
				{
					for (size_t w = 0; w < EvalPoints.getNcols(); w++)
					{
						EvalPoints.set(v, w, log(EvalPoints.get(v,w)));
					}
				}
			}
		}

		// init the basis functions' coefficient vector
		DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

		std::cout << "Grid has " << level << " Levels" << std::endl;
		std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
		std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

		// Init the grid with on payoff function
		myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);

		// Gridpoints @Money
		std::cout << "Gridpoints @Money: " << myBSSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;

		// Print the payoff function into a gnuplot file
		if (dim < 3)
		{
			myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");
			myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
			myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
			// Write interpolation error into a file
			myBSSolver->printPayoffInterpolationError2D(*alpha, "payoff_interpolation_error.grid.gnuplot", 10000, dStrike);
			if (isLogSolve == true)
			{
				myBSSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.gnuplot", true);
				myBSSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.gnuplot", false);
			}
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
			myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
		}
		else if (Solver == "AdBas")
		{
			myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
		}
		else if (Solver == "VaTim")
		{
			myBSSolver->solveVarTimestep(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
		}
		else
		{
			std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
		}

		if (dim < 3)
		{
			// Print the solved Black Scholes Equation into a gnuplot file
			myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
			myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
			myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
			if (isLogSolve == true)
			{
				myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_surplus_cart.grid.gnuplot", true);
				myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_nodal_cart.grid.gnuplot", false);
			}
		}

		// Test call @ the money
		std::vector<double> point;
		for (size_t j = 0; j < d; j++)
		{
			if (isLogSolve == true)
			{
				point.push_back(log(1.0));
			}
			else
			{
				point.push_back(1.0);
			}
		}
		std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;

		// Evaluate Cuboid
		DataVector Prices(EvalPoints.getNrows());
		myBSSolver->evaluateCuboid(*alpha, Prices, EvalPoints);
		results.push_back(Prices);

		// write solution in a additional file
		std::stringstream level_string;
		level_string << i;
		writeDataVector(Prices, tFileEvalCuboidValues+".level_"+ level_string.str());
		writeDataVector(Prices, tFileEvalCuboidValues);

		if (i > start_l)
		{
			std::cout << "=====================================================================" << std::endl;
			std::cout << "=====================================================================" << std::endl << std::endl;
			std::cout << "Calculating norms of relative errors to a grid" << std::endl;
			std::cout << "with " << i << " levels and testing-coboid" << std::endl;
			std::cout << "with the center:" << std::endl;
			for (size_t j = 0; j < d; j++)
			{
				std::cout << center[j] << " ";
			}
			std::cout << std::endl << "and " << points << " test-points in a range of " << std::endl;
			std::cout << cuboidSize*200.0 << "% per dimension:" << std::endl << std::endl;

			// Calculate relative errors and some norms
			for (size_t j = 0; j < i-start_l; j++)
			{
				DataVector maxLevel(results[i-start_l]);
				DataVector relError(results[j]);
				double maxNorm = 0.0;
				double twoNorm = 0.0;

				// calculate relative error
				relError.sub(maxLevel);
				relError.componentwise_div(maxLevel);

				// calculate max. norm of relative error
				maxNorm = relError.maxNorm();

				// calculate two norm of relative error
				twoNorm = relError.twoNorm();

				// Printing norms
				std::cout << "Level " << j + start_l << ": max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << twoNorm << std::endl;
			}
		}
		std::cout << std::endl << std::endl;

		myBSSolver->deleteGrid();
		delete alpha;

		std::cout << std::endl;
	}

	delete myEvalCuboidGen;
	delete myBSSolver;
	delete myBoundingBox;
}

/**
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call option, with Intial
 * Grid Refinement
 *
 * @param d the number of dimensions used in the Sparse Grid
 * @param l the number of levels used in the Sparse Grid
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param dStrike the strike of the option
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 * @param nIterAdaptSteps number of the iterative Grid Refinement that should be executed
 * @param dInitialAdpatDist initial distance from @the money. Is devided in every iteration by the number of the iteration
 * @param useCoarsen specifies if the grid should be coarsened between timesteps
 * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
 * @param numExecCoarsen denotes the number of coarsenings per timestep
 * @param coarsenPercent Number of removable grid points that should be tested for deletion
 * @param isLogSolve set this to true if the log-transformed Black Scholes Equation should be solved
 */
void testNUnderlyingsAdapt(size_t d, size_t l, std::string fileStoch, std::string fileBound, double dStrike, std::string payoffType,
		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, size_t nIterAdaptSteps,
		double dInitialAdpatDist, bool useCoarsen, double coarsenPercent, size_t numExecCoarsen, double coarsenThreshold, bool isLogSolve)
{
	size_t dim = d;
	size_t level = l;
	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataMatrix rho(dim,dim);

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

	// Enable Coarsening
	if (useCoarsen == true)
	{
		myBSSolver->setEnableCoarseningData(coarsenThreshold, coarsenPercent, numExecCoarsen);
	}

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
	std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

//	// refine the grid to approximate the singularity in the start solution better
//	for (size_t i = 0 ; i < nIterAdaptSteps; i++)
//	{
//		std::cout << "Refining Grid..." << std::endl;
//		myBSSolver->refineInitialGridWithPayoff(*alpha, dStrike, payoffType, (dInitialAdpatDist/(static_cast<double>(i+1))));
//		std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
//		std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
//	}

	// Generate Full Grid at @Money
	size_t oldGridSize = 0;
	size_t newGridSize = myBSSolver->getNumberGridPoints();
	size_t addedGridPoint = 0;

	size_t i = 0;
	do
	{
		oldGridSize = newGridSize;
		std::cout << "Refining Grid..." << std::endl;
		myBSSolver->refineInitialGridWithPayoffToMaxLevel(*alpha, dStrike, payoffType,  (dInitialAdpatDist/(static_cast<double>(i+1))), level);
		std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
		std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
		newGridSize = myBSSolver->getNumberGridPoints();
		addedGridPoint = newGridSize - oldGridSize;
		i++;
	} while (addedGridPoint > 0);
//
//	// Refine @Money by two Levels
//	for (size_t r = 0; r < 1; r++)
//	{
//		oldGridSize = 0;
//		newGridSize = myBSSolver->getNumberGridPoints();
//		addedGridPoint = 0;
//		do
//		{
//			oldGridSize = newGridSize;
//			std::cout << "Refining Grid..." << std::endl;
//			myBSSolver->refineInitialGridWithPayoffToMaxLevel(*alpha, dStrike, payoffType, 2.0*(1.0/((double)(1<<(level+r)))), level+1+r);
//			std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
//			std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
//			newGridSize = myBSSolver->getNumberGridPoints();
//			addedGridPoint = newGridSize - oldGridSize;
//		} while (addedGridPoint > 0);
//	}

	std::cout << std::endl << std::endl << std::endl;

	// Print the payoff function into a gnuplot file
	if (dim < 3)
	{
		myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");
		myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
		myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
		// Write interpolation error into a file
		myBSSolver->printPayoffInterpolationError2D(*alpha, "payoff_interpolation_error.grid.gnuplot", 10000, dStrike);
		if (isLogSolve == true)
		{
			myBSSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.gnuplot", true);
			myBSSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.gnuplot", false);
		}
	}

	// Gridpoints @Money
	std::cout << "Gridpoints @Money: " << myBSSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;

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
		myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
	}
	else if (Solver == "AdBas")
	{
		myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else if (Solver == "VaTim")
	{
		myBSSolver->solveVarTimestep(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else
	{
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (dim < 3)
	{
		// Print the solved Black Scholes Equation into a gnuplot file
		myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
		myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
		myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
		if (isLogSolve == true)
		{
			myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_surplus_cart.grid.gnuplot", true);
			myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_nodal_cart.grid.gnuplot", false);
		}
	}

	// Test call @ the money
	std::vector<double> point;
	for (size_t i = 0; i < d; i++)
	{
		if (isLogSolve == true)
		{
			point.push_back(log(1.0));
		}
		else
		{
			point.push_back(1.0);
		}
	}
	std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;

	// calculate relative errors
	////////////////////////////

	// read Evaluation cuboid
	DataMatrix EvalCuboid(1, dim);
	int retCuboid = readEvalutionCuboid(EvalCuboid, tFileEvalCuboid, dim);

	// read reference values for evaluation cuboid
	DataVector EvalCuboidValues(1);
	int retCuboidValues = readOptionsValues(EvalCuboidValues, tFileEvalCuboidValues);

	if (retCuboid == 0 && retCuboidValues == 0)
	{
		// If the log-transformed Black Scholes Equation is used -> transform Eval-cuboid
		if (isLogSolve == true)
		{
			for (size_t v = 0; v < EvalCuboid.getNrows(); v++)
			{
				for (size_t w = 0; w < EvalCuboid.getNcols(); w++)
				{
					EvalCuboid.set(v, w, log(EvalCuboid.get(v,w)));
				}
			}
		}

		std::cout << "Calculating relative errors..." << std::endl;
		// Evaluate Cuboid
		DataVector Prices(EvalCuboid.getNrows());
		myBSSolver->evaluateCuboid(*alpha, Prices, EvalCuboid);

		DataVector relError(Prices);
		double maxNorm = 0.0;
		double twoNorm = 0.0;

		// calculate relative error
		relError.sub(EvalCuboidValues);
		relError.componentwise_div(EvalCuboidValues);

		// calculate max. norm of relative error
		maxNorm = relError.maxNorm();

		// calculate two norm of relative error
		twoNorm = relError.twoNorm();

		// Printing norms
		std::cout << "Results: max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << twoNorm << std::endl;
	}
	else
	{
		std::cout << "Couldn't open evaluation cuboid data -> skipping tests!" << std::endl << std::endl;
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
 * @param dStrike the strike of the option
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 * @param refinePercent percantage of points that should be refined before Black Scholes Equation is solved
 * @param nIterAdaptSteps number of the iterative Grid Refinement that should be executed
 * @param dRefineThreshold Threshold for a point's surplus for refining this point
 * @param useCoarsen specifies if the grid should be coarsened between timesteps
 * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
 * @param numExecCoarsen denotes the number of coarsenings per timestep
 * @param coarsenPercent Number of removable grid points that should be tested for deletion
 * @param nIterAdaptSteps number of the iterative Grid Refinement that should be executed
 * @param isLogSolve set this to true if the log-transformed Black Scholes Equation should be solved
 */
void testNUnderlyingsAdaptSurplus(size_t d, size_t l, std::string fileStoch, std::string fileBound, double dStrike,
		std::string payoffType, double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps,
		std::string Solver, double refinePercent, size_t nIterAdaptSteps, double dRefineThreshold,
		bool useCoarsen, double coarsenPercent, size_t numExecCoarsen, double coarsenThreshold, bool isLogSolve)
{
	size_t dim = d;
	size_t level = l;
	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;

	DataVector mu(dim);
	DataVector sigma(dim);
	DataMatrix rho(dim,dim);

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

	// Enable Coarsening
	if (useCoarsen == true)
	{
		myBSSolver->setEnableCoarseningData(coarsenThreshold, coarsenPercent, numExecCoarsen);
	}

	// init the basis functions' coefficient vector
	DataVector* alpha = new DataVector(myBSSolver->getNumberGridPoints());

	std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
	std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

	// Init the grid with on payoff function
	myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);

	// refine the grid to approximate the singularity in the start solution better
	for (size_t i = 0 ; i < nIterAdaptSteps; i++)
	{
		std::cout << "Refining Grid..." << std::endl;
		myBSSolver->refineInitialGridSurplus(*alpha, refinePercent, dRefineThreshold);
		myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
		std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
		std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
	}
	std::cout << std::endl << std::endl << std::endl;

//	size_t oldGridSize = 0;
//	size_t newGridSize = myBSSolver->getNumberGridPoints();
//	size_t addedGridPoint = 0;
//	do
//	{
//		oldGridSize = newGridSize;
//		std::cout << "Refining Grid..." << std::endl;
//		myBSSolver->refineSurplusToMaxLevel(*alpha, dRefineThreshold, level);
//		myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//		std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
//		std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
//		newGridSize = myBSSolver->getNumberGridPoints();
//		addedGridPoint = newGridSize - oldGridSize;
//	} while (addedGridPoint > 0);
//	std::cout << std::endl << std::endl << std::endl;

//	for (size_t i = 0 ; i < nIterAdaptSteps; i++)
//	{
//		std::cout << "Refining Grid..." << std::endl;
//		myBSSolver->refineSurplusToMaxLevel(*alpha, dRefineThreshold, level);
//		myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//		std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
//		std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
//	}
//	std::cout << std::endl << std::endl << std::endl;

	// Print the payoff function into a gnuplot file
	if (dim < 3)
	{
		myBSSolver->printGrid(*alpha, 20, "payoff.gnuplot");
	}
	myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
	myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
	// Write interpolation error into a file
	myBSSolver->printPayoffInterpolationError2D(*alpha, "payoff_interpolation_error.grid.gnuplot", 10000, dStrike);
	if (isLogSolve == true)
	{
		myBSSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.gnuplot", true);
		myBSSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.gnuplot", false);
	}

	// Gridpoints @Money
	std::cout << "Gridpoints @Money: " << myBSSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;

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
		myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
	}
	else if (Solver == "AdBas")
	{
		myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else if (Solver == "VaTim")
	{
		myBSSolver->solveVarTimestep(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else
	{
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (dim < 3)
	{
		// Print the solved Black Scholes Equation into a gnuplot file
		myBSSolver->printGrid(*alpha, 20, "solvedBS.gnuplot");
	}
	myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
	myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
	if (isLogSolve == true)
	{
		myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_surplus_cart.grid.gnuplot", true);
		myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_nodal_cart.grid.gnuplot", false);
	}

	std::vector<double> point;
	for (size_t i = 0; i < d; i++)
	{
		if (isLogSolve == true)
		{
			point.push_back(log(1.0));
		}
		else
		{
			point.push_back(1.0);
		}
	}
	std::cout << "Optionprice at testpoint: " << myBSSolver->evaluatePoint(point, *alpha) << std::endl << std::endl;

	// calculate relative errors
	////////////////////////////

	// read Evaluation cuboid
	DataMatrix EvalCuboid(1, dim);
	int retCuboid = readEvalutionCuboid(EvalCuboid, tFileEvalCuboid, dim);

	// read reference values for evaluation cuboid
	DataVector EvalCuboidValues(1);
	int retCuboidValues = readOptionsValues(EvalCuboidValues, tFileEvalCuboidValues);

	if (retCuboid == 0 && retCuboidValues == 0)
	{
		// If the log-transformed Black Scholes Equation is used -> transform Eval-cuboid
		if (isLogSolve == true)
		{
			for (size_t v = 0; v < EvalCuboid.getNrows(); v++)
			{
				for (size_t w = 0; w < EvalCuboid.getNcols(); w++)
				{
					EvalCuboid.set(v, w, log(EvalCuboid.get(v,w)));
				}
			}
		}

		std::cout << "Calculating relative errors..." << std::endl;
		// Evaluate Cuboid
		DataVector Prices(EvalCuboid.getNrows());
		myBSSolver->evaluateCuboid(*alpha, Prices, EvalCuboid);

		DataVector relError(Prices);
		double maxNorm = 0.0;
		double twoNorm = 0.0;

		// calculate relative error
		relError.sub(EvalCuboidValues);
		relError.componentwise_div(EvalCuboidValues);

		// calculate max. norm of relative error
		maxNorm = relError.maxNorm();

		// calculate two norm of relative error
		twoNorm = relError.twoNorm();

		// Printing norms
		std::cout << "Results: max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << twoNorm << std::endl;
	}
	else
	{
		std::cout << "Couldn't open evaluation cuboid data -> skipping tests!" << std::endl << std::endl;
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
	DataMatrix rho(dim, dim);
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
		myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
	}
	else if (Solver == "AdBas")
	{
		myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else if (Solver == "VaTim")
	{
		myBSSolver->solveVarTimestep(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else
	{
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic" || Solver == "AdBas" || Solver == "VaTim")
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
	mySStream << "  solveND             Solves an European Call/Put option" << std::endl;
	mySStream << "                      for N assets on a regular sparse grid" << std::endl << std::endl;
	mySStream << "  solveNDanalyze      same as solveND, but the option is" << std::endl;
	mySStream << "                      solved for several regular grids with" << std::endl;
	mySStream << "                      different numbers of levels" << std::endl << std::endl;
	mySStream << "  solveNDadapt        Solves an European Call/Put option" << std::endl;
	mySStream << "                      on an initial refined grid (based on a" << std::endl;
	mySStream << "                      regular grid) by analyzing the payoff" << std::endl << std::endl;
	mySStream << "  solveNDadaptSurplus Solves an European Call/Up option" << std::endl;
	mySStream << "                      on an initial refined grid based on" << std::endl;
	mySStream << "                      the hierarchical surplus" << std::endl << std::endl;
	mySStream << "  solveBonn  Solves an option delivered in Bonn's format" << std::endl << std::endl << std::endl;

	mySStream << "Several files are needed to specify inputs:" << std::endl;
	mySStream << "-----------------------------------------------------" << std::endl;
	mySStream << "file_Boundaries:  this file contains the grid's bounding box" << std::endl;
	mySStream << "                  for every dimension this file contains a" << std::endl;
	mySStream << "                  tuple with the boundaries." << std::endl;
	mySStream << "Example (3 dimensions):" << std::endl;
	mySStream << "                  0.0 2.5" << std::endl;
	mySStream << "                  0.0 2.5" << std::endl;
	mySStream << "                  0.0 2.5" << std::endl << std::endl << std::endl;

	mySStream << "file_Stochdata:   this file contains the option's asset's" << std::endl;
	mySStream << "                  expected values, standard deviations and" << std::endl;
	mySStream << "                  correlations. The i-th line contains" << std::endl;
	mySStream << "                  followning data:" << std::endl;
	mySStream << "                  mu_i sigma_i rho_i_0 ... rhi_i_d" << std::endl;
	mySStream << "Example (3 dimensions):" << std::endl;
	mySStream << "                  0.05 0.4 1.0 0.1 0.2" << std::endl;
	mySStream << "                  0.05 0.5 0.1 1.0 0.3" << std::endl;
	mySStream << "                  0.05 0.6 0.2 0.3 1.0" << std::endl << std::endl << std::endl;

	mySStream << "file_analyze:     this file contains the options for" << std::endl;
	mySStream << "                  the analyzing runs. This file contains" << std::endl;
	mySStream << "                  two lines: The first lines is the center" << std::endl;
	mySStream << "                  of the evaluation cuboid. The second one" << std::endl;
	mySStream << "                  contains two values: The first one is the" << std::endl;
	mySStream << "                  the evaluations cuboid's expansion in percent" << std::endl;
	mySStream << "                  around the former specified center. The " << std::endl;
	mySStream << "                  second one is the number of points" << std::endl;
	mySStream << "                  in every dimension in the evaluation" << std::endl;
	mySStream << "                  cuboid." << std::endl;
	mySStream << "Example (3 dimensions):" << std::endl;
	mySStream << "                  1.0 1.0 1.0" << std::endl;
	mySStream << "                  0.05 20" << std::endl << std::endl << std::endl;

	mySStream << "Execution modes descriptions:" << std::endl;
	mySStream << "-----------------------------------------------------" << std::endl;
	mySStream << "solveND" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	Strike: the strike" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 cart" << std::endl;
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
	mySStream << "	Strikes: the strike" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	file_analyze: file containing the analyzing options" << std::endl;
	mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 2 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 anal.data cart" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
	mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
	mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
	mySStream << "For all cases following files are generated:" << std::endl;
	mySStream << "	EvalCuboidPoints.data: containing the evaluation" << std::endl;
	mySStream << "		cuboid" << std::endl;
	mySStream << "	HighLevelOptionValue.data: containing the option's" << std::endl;
	mySStream << "		for the highest leveled grid." << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "solveNDadapt" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	Strike: the strike" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
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
	mySStream << "3 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 5 0.5" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
	mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
	mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "solveNDadaptFull" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	Strike: the strike" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
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
	mySStream << "	Coarsening Percent: Percent of gird point that should be removed" << std::endl;
	mySStream << "	Num Coarsenings: Number of coarsenings per timestep" << std::endl;
	mySStream << "	Coarsening Threshold: Threshold of point's surplus to remove point" << std::endl;
	mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 5 0.5 80 2 1e-6 cart" << std::endl;
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
	mySStream << "	Strike: the strike" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	Adapt-Refinement Percent: Percent of grid points that should be refined" << std::endl;
	mySStream << "	numAdaptRefinement: Number of adaptive refinements at the beginning" << std::endl;
	mySStream << "	refinementThreshold: Threshold of point's surplus to refine point" << std::endl;
	mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 50 10 1e-10 cart" << std::endl;
	mySStream << std::endl;
	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	mySStream << "And for dim>2 Bonn formated Sparse Grid files:" << std::endl;
	mySStream << "	payoff_Nd.bonn: the start condition" << std::endl;
	mySStream << "	solvedBS_Nd.bonn: the numerical solution" << std::endl;
	mySStream << std::endl << std::endl;

	mySStream << "solveNDadaptSurplusFull" << std::endl << "------" << std::endl;
	mySStream << "the following options must be specified:" << std::endl;
	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	mySStream << "	Strike: the strike" << std::endl;
	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	mySStream << "	r: the riskfree rate" << std::endl;
	mySStream << "	T: time to maturity" << std::endl;
	mySStream << "	dT: timestep size" << std::endl;
	mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	mySStream << "	Adapt-Refinement Percent: Percent of grid points that should be refined" << std::endl;
	mySStream << "	numAdaptRefinement: Number of adaptive refinements at the beginning" << std::endl;
	mySStream << "	refinementThreshold: Threshold of point's surplus to refine point" << std::endl;
	mySStream << "	Coarsening Percent: Percent of gird point that should be removed" << std::endl;
	mySStream << "	Num Coarsenings: Number of coarsenings per timestep" << std::endl;
	mySStream << "	Coarsening Threshold: Threshold of point's surplus to remove point" << std::endl;
	mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords" << std::endl;
	mySStream << std::endl;
	mySStream << "Example:" << std::endl;
	mySStream << "3 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 50 10 1e-10 80 2 1e-6 cart" << std::endl;
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

	if (option == "solveND")
	{
		if (argc != 15)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			double dStrike;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[5]);
			fileBound.assign(argv[4]);
			dStrike = atof(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[11]);

			std::string coordsType;
			bool coords = false;
			coordsType.assign(argv[14]);
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

			testNUnderlyings(atoi(argv[2]), atoi(argv[3]), fileStoch, fileBound, dStrike, payoff, atof(argv[8]), (size_t)(atof(argv[9])/atof(argv[10])), atof(argv[10]), atoi(argv[12]), atof(argv[13]), solver, coords);
		}
	}
	else if (option == "solveNDanalyze")
	{
		if (argc != 17)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			double dStrike;
			std::string fileAnalyze;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[6]);
			fileBound.assign(argv[5]);
			dStrike = atof(argv[7]);
			fileAnalyze.assign(argv[15]);
			payoff.assign(argv[8]);
			solver.assign(argv[12]);

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


			testNUnderlyingsAnalyze(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), fileStoch, fileBound, dStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, fileAnalyze, coords);
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
			double dStrike;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[5]);
			fileBound.assign(argv[4]);
			dStrike = atof(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[11]);

//			std::string coordsType;
//			bool coords = false;
//			coordsType.assign(argv[16]);
//			if (coordsType == "cart")
//			{
//				coords = false;
//			}
//			else if (coordsType == "log")
//			{
//				coords = true;
//			}
//			else
//			{
//				std::cout << "Unsupported coordinate option! cart or log are supported!" << std::endl;
//				std::cout << std::endl << std::endl;
//				writeHelp();
//			}

			testNUnderlyingsAdapt(atoi(argv[2]), atoi(argv[3]), fileStoch, fileBound, dStrike, payoff, atof(argv[8]), (size_t)(atof(argv[9])/atof(argv[10])), atof(argv[10]), atoi(argv[12]), atof(argv[13]), solver, atoi(argv[14]), atof(argv[15]), false, 0.0, 0, 0.0, false);
		}
	}
	else if (option == "solveNDadaptSurplus")
	{
		if (argc != 18)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			double dStrike;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[5]);
			fileBound.assign(argv[4]);
			dStrike = atof(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[11]);

			std::string coordsType;
			bool coords = false;
			coordsType.assign(argv[17]);
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

			testNUnderlyingsAdaptSurplus(atoi(argv[2]), atoi(argv[3]), fileStoch, fileBound, dStrike, payoff, atof(argv[8]), (size_t)(atof(argv[9])/atof(argv[10])), atof(argv[10]), atoi(argv[12]), atof(argv[13]), solver, atoi(argv[14]), atoi(argv[15]), atof(argv[16]), false, 0.0, 0, 0.0, coords);
		}
	}
	else if (option == "solveNDadaptFull")
	{
		if (argc != 20)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			double dStrike;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[5]);
			fileBound.assign(argv[4]);
			dStrike = atof(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[11]);

			std::string coordsType;
			bool coords = false;
			coordsType.assign(argv[19]);
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

			testNUnderlyingsAdapt(atoi(argv[2]), atoi(argv[3]), fileStoch, fileBound, dStrike, payoff, atof(argv[8]), (size_t)(atof(argv[9])/atof(argv[10])), atof(argv[10]), atoi(argv[12]), atof(argv[13]), solver, atoi(argv[14]), atof(argv[15]), true, atof(argv[16]), atoi(argv[17]), atof(argv[18]), coords);
		}
	}
	else if (option == "solveNDadaptSurplusFull")
	{
		if (argc != 21)
		{
			writeHelp();
		}
		else
		{
			std::string fileStoch;
			std::string fileBound;
			double dStrike;
			std::string ani;
			std::string solver;
			std::string payoff;

			fileStoch.assign(argv[5]);
			fileBound.assign(argv[4]);
			dStrike = atof(argv[6]);
			payoff.assign(argv[7]);
			solver.assign(argv[11]);

			std::string coordsType;
			bool coords = false;
			coordsType.assign(argv[20]);
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

			testNUnderlyingsAdaptSurplus(atoi(argv[2]), atoi(argv[3]), fileStoch, fileBound, dStrike, payoff, atof(argv[8]), (size_t)(atof(argv[9])/atof(argv[10])), atof(argv[10]), atoi(argv[12]), atof(argv[13]), solver, atoi(argv[14]), atoi(argv[15]), atof(argv[16]), true, atof(argv[17]), atoi(argv[18]), atof(argv[19]), coords);
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
