/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "sgpp_base.hpp"
#include "sgpp_pde.hpp"
#include "sgpp_finance.hpp"
#include "sgpp_parallel.hpp"
#include "sgpp_solver.hpp"
#include "sgpp_datadriven.hpp"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <limits>

// default number of Implicit Euler steps before starting with Crank Nicolson approach
#define CRNIC_IMEUL_STEPS 3
// default value for epsilon in gridpoints @money
#define DFLT_EPS_AT_MONEY 0.0
// Resolution (num points per dimension) of the created gnuplots
#define PLOT_RESOLUTION 40

using namespace sg;
using namespace sg::pde;
using namespace std;

double alphaDone;
double vProbe;
double sProbe;
double refinementThresh;
size_t numGridPoints;

/**
 * Calls the writeHelp method in the HestonSolver Object
 * after creating a screen.
 */
void writeHelp()
{
	//	sg::finance::HestonSolver* myHestonSolver = new sg::finance::HestonSolver();
	//
	//	myHestonSolver->initScreen();
	//
	//	delete myHestonSolver;
	//
	//	std::stringstream mySStream;
	//
	//	mySStream << "Some instructions for the use of Black Scholes Solver:" << std::endl;
	//	mySStream << "------------------------------------------------------" << std::endl << std::endl;
	//	mySStream << "Available execution modes are:" << std::endl;
	//	mySStream << "  solveND             Solves an European Call/Put option" << std::endl;
	//	mySStream << "                      for N assets on a regular sparse grid" << std::endl << std::endl;
	//	mySStream << "  solveNDanalyze      same as solveND, but the option is" << std::endl;
	//	mySStream << "                      solved for several regular grids with" << std::endl;
	//	mySStream << "                      different numbers of levels" << std::endl << std::endl;
	//	mySStream << "  solveNDadaptSurplus Solves an European Call/Up option" << std::endl;
	//	mySStream << "                      on a refined grid based on" << std::endl;
	//	mySStream << "                      the hierarchical surplus" << std::endl << std::endl;
	//	mySStream << "  solveNDadaptSurplusSubDomain   Same as above but" << std::endl;
	//	mySStream << "						a normal distribution is used" << std::endl;
	//	mySStream << "						to do refinement just near the strike!" << std::endl << std::endl;
	//	mySStream << "  solveBonn  Solves an option delivered in Bonn's format" << std::endl << std::endl << std::endl;
	//
	//	mySStream << "Several files are needed to specify inputs:" << std::endl;
	//	mySStream << "-----------------------------------------------------" << std::endl;
	//	mySStream << "file_Boundaries:  this file contains the grid's bounding box" << std::endl;
	//	mySStream << "                  for every dimension this file contains a" << std::endl;
	//	mySStream << "                  tuple with the boundaries." << std::endl;
	//	mySStream << "Example (3 dimensions):" << std::endl;
	//	mySStream << "                  0.0 2.5" << std::endl;
	//	mySStream << "                  0.0 2.5" << std::endl;
	//	mySStream << "                  0.0 2.5" << std::endl << std::endl << std::endl;
	//
	//	mySStream << "file_Stochdata:   this file contains the option's asset's" << std::endl;
	//	mySStream << "                  expected values, standard deviations and" << std::endl;
	//	mySStream << "                  correlations. The i-th line contains" << std::endl;
	//	mySStream << "                  followning data:" << std::endl;
	//	mySStream << "                  mu_i sigma_i rho_i_0 ... rhi_i_d" << std::endl;
	//	mySStream << "Example (3 dimensions):" << std::endl;
	//	mySStream << "                  0.05 0.4 1.0 0.1 0.2" << std::endl;
	//	mySStream << "                  0.05 0.5 0.1 1.0 0.3" << std::endl;
	//	mySStream << "                  0.05 0.6 0.2 0.3 1.0" << std::endl << std::endl << std::endl;
	//
	//	mySStream << "file_analyze:     this file contains the options for" << std::endl;
	//	mySStream << "                  the analyzing runs. This file contains" << std::endl;
	//	mySStream << "                  two parts: The first lines is the " << std::endl;
	//	mySStream << "                  evaluation cuboid as bounding box. " << std::endl;
	//	mySStream << "                  The second one is the number of points" << std::endl;
	//	mySStream << "                  in every dimension in the evaluation" << std::endl;
	//	mySStream << "                  cuboid." << std::endl;
	//	mySStream << "Example (3 dimensions):" << std::endl;
	//	mySStream << "                  0.0 1.0" << std::endl;
	//	mySStream << "                  0.0 1.0" << std::endl;
	//	mySStream << "                  0.0 1.0" << std::endl;
	//	mySStream << "                  20" << std::endl << std::endl << std::endl;
	//
	//	mySStream << "Execution modes descriptions:" << std::endl;
	//	mySStream << "-----------------------------------------------------" << std::endl;
	//	mySStream << "solveND" << std::endl << "------" << std::endl;
	//	mySStream << "the following options must be specified:" << std::endl;
	//	mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords; PAT principal axis transform" << std::endl;
	//	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	//	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	//	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	//	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	//	mySStream << "	Strike: the strike" << std::endl;
	//	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	//	mySStream << "	r: the riskfree rate" << std::endl;
	//	mySStream << "	T: time to maturity" << std::endl;
	//	mySStream << "	dT: timestep size" << std::endl;
	//	mySStream << "	Solver: the solver to use: ExEul, ImEul, CrNic, AdBas, SC2:epsilon:gamma, SCH:epsilon:gamma or SCI:epsilon:sc:gamma" << std::endl;
	//	mySStream << "          (for explanations of the options, see end of help!)" << std::endl;
	//	mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
	//	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	//	mySStream << std::endl;
	//	mySStream << "Example:" << std::endl;
	//	mySStream << "cart 3 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001" << std::endl;
	//	mySStream << std::endl;
	//	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	//	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	//	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	//	mySStream << std::endl << std::endl;
	//
	//	mySStream << "solveNDanalyze" << std::endl << "------" << std::endl;
	//	mySStream << "the following options must be specified:" << std::endl;
	//	mySStream << "	Coordinates: cart: cartesian coordinates; log: log coords; PAT principal axis transform" << std::endl;
	//	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	//	mySStream << "	level_start: number of levels within the Sparse Grid (start)" << std::endl;
	//	mySStream << "	level_end: number of levels within the Sparse Grid (end)" << std::endl;
	//	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	//	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	//	mySStream << "	Strikes: the strike" << std::endl;
	//	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	//	mySStream << "	r: the riskfree rate" << std::endl;
	//	mySStream << "	T: time to maturity" << std::endl;
	//	mySStream << "	dT: timestep size" << std::endl;
	//	mySStream << "	Solver: the solver to use: ExEul, ImEul, CrNic, AdBas, SC2:epsilon:gamma, SCH:epsilon:gamma or SCI:epsilon:sc:gamma" << std::endl;
	//	mySStream << "          (for explanations of the options, see end of help!)" << std::endl;
	//	mySStream << "	CGIterations: Maxmimum number of iterations used in CG method" << std::endl;
	//	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	//	mySStream << "	file_analyze: file containing the analyzing options" << std::endl;
	//	mySStream << std::endl;
	//	mySStream << "Example:" << std::endl;
	//	mySStream << "cart 3 2 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 anal.data" << std::endl;
	//	mySStream << std::endl;
	//	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	//	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	//	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	//	mySStream << "For all cases following files are generated:" << std::endl;
	//	mySStream << "	EvalCuboidPoints.data: containing the evaluation" << std::endl;
	//	mySStream << "		cuboid" << std::endl;
	//	mySStream << "	HighLevelOptionValue.data: containing the option's" << std::endl;
	//	mySStream << "		for the highest leveled grid." << std::endl;
	//	mySStream << std::endl << std::endl;
	//
	//	mySStream << "solve1Danalyze" << std::endl << "------" << std::endl;
	//	mySStream << "compares with analytic solution" << std::endl;
	//	mySStream << "same parameter configuration as in mode 'solveNDanalyze', but dim has to be set to 1 and mu=r" << std::endl;
	//	mySStream << std::endl << std::endl;
	//
	//	mySStream << "solveNDadaptSurplus/solveNDadaptSurplusSubDomain" << std::endl << "------" << std::endl;
	//	mySStream << "the following options must be specified:" << std::endl;
	//	mySStream << "	Coordinates: cart: cartesian coordinates; log: log coords; PAT principal axis transform" << std::endl;
	//	mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
	//	mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
	//	mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
	//	mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
	//	mySStream << "	Strike: the strike" << std::endl;
	//	mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
	//	mySStream << "	r: the riskfree rate" << std::endl;
	//	mySStream << "	T: time to maturity" << std::endl;
	//	mySStream << "	dT: timestep size" << std::endl;
	//	mySStream << "	Solver: the solver to use: ExEul, ImEul, CrNic, AdBas, SC2:epsilon:gamma, SCH:epsilon:gamma or SCI:epsilon:sc:gamma" << std::endl;
	//	mySStream << "          (for explanations of the options, see end of help!)" << std::endl;
	//	mySStream << "	CGIterations: Maxmimum number of iterations used in CG method" << std::endl;
	//	mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
	//	mySStream << "	RefinementMode: classic or maxLevel" << std::endl;
	//	mySStream << "	MaxRefinement Level: Max. Level for refinement" << std::endl;
	//	mySStream << "	numAdaptRefinement: Number of adaptive refinements at the beginning" << std::endl;
	//	mySStream << "	refinementThreshold: Threshold of point's surplus to refine point" << std::endl;
	//	mySStream << "	adapt-mode during solving: none, coarsen, refine, coarsenNrefine" << std::endl;
	//	mySStream << "	Coarsening Threshold: Threshold of point's surplus to remove point" << std::endl;
	//	mySStream << std::endl;
	//	mySStream << "Example:" << std::endl;
	//	mySStream << "cart 3 5 " << "bound.data stoch.data 1.0 std_euro_call "<< "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 classic 0 5 1e-10 coarsen 1e-6" << std::endl;
	//	mySStream << std::endl;
	//	mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
	//	mySStream << "	payoff.gnuplot: the start condition" << std::endl;
	//	mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
	//	mySStream << std::endl << std::endl;
	//
	//	mySStream << "analyzeTimeStepping" << std::endl << "------" << std::endl;
	//	mySStream << "for analyzing time-stepping of SCAC, SCH, SCBDF or SCEJ" << std::endl;
	//	mySStream << "same function call as in 'solveNDanalyze' plus additional parameters: " << std::endl;
	//	mySStream << "epsTime: epsilon for timestepping " << std::endl;
	//	mySStream << "cTime: constant for timestepping (optional, only in case of SCEJ needed)" << std::endl;
	//	mySStream << std::endl << std::endl;
	//
	//	mySStream << "options for time-stepping:" << std::endl << "------" << std::endl;
	//	mySStream << "   * ExEul                  explicit Euler" << std::endl;
	//	mySStream << "   * ImEul                  implicit Euler" << std::endl;
	//	mySStream << "   * CrNic                  Crank-Nicoloson" << std::endl;
	//	mySStream << "   * AdBas                  Adams-Bashforth (explicit 2-step method)" << std::endl;
	//	mySStream << " adaptive time-steppings:" << std::endl;
	//	mySStream << "   * SC2:epsilon:gamma      predictor-corrector scheme (Crank-Nicolson/Adams-Bashforth)" << std::endl;
	//	mySStream << "   * SCH:epsilon:gamma      step-doubling (Crank-Nicolson)" << std::endl;
	//	mySStream << "   * SCI:epsilon:sc:gamma   increase-in-value (Crank-Nicolson)" << std::endl;
	//	mySStream << " hereby: " << std::endl;
	//	mySStream << "   epsilon = rel. discretization accuracy " << std::endl;
	//	mySStream << "             (typically eps = 5e-3 to 1e-4)" << std::endl;
	//	mySStream << "   gamma   = factor by which the step size is decreased (dt = gamma*dt) in case of rejection" << std::endl;
	//	mySStream << "             (typically gamma = 0.5 to 0.9)" << std::endl;
	//	mySStream << "   sc      = specific factor of increase-in-value for limiting the increase" << std::endl;
	//	mySStream << "             (typically sc = 1 or 10)" << std::endl;
	//	mySStream << std::endl << std::endl;
	//
	//	mySStream << std::endl << std::endl;
	//	std::cout << mySStream.str() << std::endl;
}

/**
 * reads the values of mu, sigma and rho of all assets from
 * a file and stores them into three separated DataVectors
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssets the of Assets stored in the file
 * @param mu DataVector for the exspected values
 * @param sigma DataVector for standard deviation
 * @param rho DataMatrix for the correlations
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readStochasticData(std::string tFile, size_t numAssets, sg::base::DataVector& xi, sg::base::DataVector& theta, sg::base::DataVector& kappa, sg::base::DataMatrix& hMatrix)
{
	// For the Heston model we need the following stochastic process data
	// mu: the drift vector, one value for each asset
	// v: the variance vector, one value for each asset
	// rho: each rho value represents the correlation between the stock prices process and the variance process for an individual asset. For d assets we have a vector of size d here.
	// xi: each xi value represents the volatility of the volatility (also called volatility of the variance) for a particular asset. For d assets we have a vector of size d here.
	// theta: each theta value represents the long-run variance for a particular asset. For d assets we have a vector of size d here.
	// kappa: each kappa value represents the mean-reversion rate (i.e. the speed which the variance goes back to its long-run value whenever it's not equal to it). For d assets we have a vector of size d here.
	// H: the correlation matrix, of size dxd, where d is the number of assets

	// Explanation: Correlation in the Heston model:
	// For one asset, we have two processes, dS and dv, the Wiener processes of which are correlated by a factor p.
	// So that means that for one asset we essentially have a 2x2 H matrix. The h11 entry of this matrix is the correlation between the S Wiener process and the
	// S Wiener process (this is usually one). The h12 entry of the matrix is the correlation between the S Wiener process and the v wiener process, and so on for h21 and h22.
	// For two assets h would be 4x4 and so on.
	// Example of the file for one and two assets:
	//
	// (for one asset)
	// xi theta kappa h11 h12 h21 h22
	//
	// (for two assets)
	// xi1 theta1 kappa1 h11 h12 h13 h14 h21 h22 h23 h24
	// xi2 theta2 kappa2 h31 h32 h33 h34 h41 h42 h43 h44
	//
	// ...and so on for more assets

	std::fstream file;
	double cur_xi;
	double cur_theta;
	double cur_kappa;
	double cur_h;
	size_t hMatrixDim = 2*numAssets;

	// Count the number of elements in the file
	file.open(tFile.c_str());
	if(!file.is_open())
	{
		std::cout << "Error cannot read file: " << tFile << std::endl;
		return -1;
	}

	// There must be a xi, theta and kappa for each asset
	// In addition, there must be a total of (2xnumAssets)^2 hmatrix elements
	// Thus, the element count must be (3xnumAssets)+(2xnumAssets)^2
	size_t t = 0;
	double test;
	do
	{
		file >> test;
		t++;
	} while (!file.eof());
	file.close();
	if (t < ((hMatrixDim*hMatrixDim)+(3*numAssets)))
	{
		std::cout << "Invalid stoch file: " << tFile << " Last Value:" << test << std::endl;
		return -1;
	}

	// Read the file (outermost loop here iterates over each asset line of text)
	file.open(tFile.c_str());
	for (size_t i = 0; i < numAssets; i++)
	{
		// Read in the xi, theta and kappa values for asset i from the file
		file >> cur_xi;
		file >> cur_theta;
		file >> cur_kappa;

		// Set the xi, theta and kappa values for asset i
		xi.set(i, cur_xi);
		theta.set(i, cur_theta);
		kappa.set(i, cur_kappa);

		// Now we deal with the h matrix values.
		// On each text line (i.e. for each asset), the number of matrix entries is equal
		// to 2*(2*numAssets), i.e. two matrix rows of (2*numAssets) each.
		// We will process each of these rows in a separate loop. The first loop handles the first row (i.e. row index 2*i)
		// and the second loop handles the second row (i.e. row index 2*i + 1).
		for (size_t j = 0; j < hMatrixDim; j++)
		{
			file >> cur_h;
			hMatrix.set(2*i,j, cur_h);
		}
		for (size_t j = 0; j < hMatrixDim; j++)
		{
			file >> cur_h;
			hMatrix.set(2*i+1,j, cur_h);
		}
	}

	file.close();

	return 0;
}

/**
 * reads the values of the Bounding Box
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssets the of Assets stored in the file
 * @param BoundaryArray Pointer to the Bounding Box array
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readBoudingBoxData(std::string tFile, size_t numDims, sg::base::DimensionBoundary* BoundaryArray)
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

	// Get number of elements in bound file, must be 2*numDims (i.e. an upper and lower bound for each dimension)
	size_t j = 0;
	double test;
	do
	{
		file >> test;
		j++;
	} while (!file.eof());
	file.close();
	if (j < (numDims*2))
	{
		std::cout << "Invalid boundary file (j=" << j << "): " << tFile << " Last Value:" << test << std::endl;
		return -1;
	}

	// Read the boundary data...two elements (left and right) for each dimension.
	file.open(tFile.c_str());
	for (size_t i = 0; i < numDims; i++)
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
 * @param coordsType set the type of coordinates that should be used: cart, log, PAT
 */
void testNUnderlyings(size_t numAssets, size_t l, std::string fileStoch, std::string fileBound, double dStrike, std::string payoffType,
		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, std::string coordsType)
{
	size_t numberOfAssets = numAssets;
	size_t pdeDim = numberOfAssets*2;
	size_t level = l;
	size_t timesteps = timeSt;
	double stepsize = dt;
	size_t CGiterations = CGIt;
	double CGepsilon = CGeps;
	double maxStock = 0.0;

	sg::base::DataVector theta(numberOfAssets);
	sg::base::DataVector xi(numberOfAssets);
	sg::base::DataVector kappa(numberOfAssets);
	sg::base::DataMatrix hMatrix(pdeDim,pdeDim);

	double r = riskfree;

	if (readStochasticData(fileStoch, numberOfAssets, xi, theta, kappa, hMatrix) != 0)
	{
		return;
	}

	// We have boundary data for each dimension in the PDE
	sg::base::DimensionBoundary* myBoundaries = new sg::base::DimensionBoundary[pdeDim];
	if (readBoudingBoxData(fileBound, pdeDim, myBoundaries) != 0)
	{
		return;
	}

	sg::finance::HestonSolver* myHestonSolver;
	if (coordsType == "log")
	{
		myHestonSolver = new sg::finance::HestonSolver(true);
	}
	else if (coordsType == "cart")
	{
		myHestonSolver = new sg::finance::HestonSolver(false);
	}
	else
	{
		// Write Error
		std::cout << "Unsupported grid transformation!" << std::endl;
		std::cout << std::endl << std::endl;
		writeHelp();
	}

	sg::base::BoundingBox* myBoundingBox = new sg::base::BoundingBox(pdeDim, myBoundaries);
	//	if (numberOfAssets == 1)
	//	{
	//		maxStock = myBoundaries[0].rightBoundary;
	//	}
	delete[] myBoundaries;

	// init Screen Object
	myHestonSolver->initScreen();

	// Construct a grid
	myHestonSolver->constructGrid(*myBoundingBox, level);

	std::string adaptSolvingMode = "refine";
	std::string refinementMode = "classic";
	size_t maxRefineLevel = 11;
	double coarsenThreshold = 0.0;
	double dRefineThreshold = 0.0000001;// See Alex's second thesis
//	double dRefineThreshold = refinementThresh;

	// Set coarsening dat
	//	myHestonSolver->setEnableCoarseningData(adaptSolvingMode, refinementMode, maxRefineLevel, -1, coarsenThreshold, dRefineThreshold);

	// init the basis functions' coefficient vector
	sg::base::DataVector* alpha = new sg::base::DataVector(myHestonSolver->getNumberGridPoints());

	std::cout << "Grid has " << level << " Levels" << std::endl;
	std::cout << "Initial Grid size: " << myHestonSolver->getNumberGridPoints() << std::endl;
	std::cout << "Initial Grid size (inner): " << myHestonSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

	size_t nIterAdaptSteps = 5;
	bool useNormalDist = true;
	int numRefinePoints = -1;
	std::vector<double> norm_mu;
	std::vector<double> norm_sigma;
	norm_mu.push_back(0.5); norm_mu.push_back(5);
	norm_sigma.push_back(0.5); norm_sigma.push_back(5);

	// refine the grid to approximate the singularity in the start solution better
//	if (refinementMode == "classic")
//	{
//		for (size_t i = 0 ; i < nIterAdaptSteps; i++)
//		{
//			std::cout << "Refining Grid..." << std::endl;
//			if (useNormalDist == true)
//			{
//				myHestonSolver->refineInitialGridSurplusSubDomain(*alpha, numRefinePoints, dRefineThreshold, norm_mu, norm_sigma);
//			}
//			else
//			{
//				myHestonSolver->refineInitialGridSurplus(*alpha, numRefinePoints, dRefineThreshold);
//			}
//			myHestonSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//			std::cout << "Refined Grid size: " << myHestonSolver->getNumberGridPoints() << std::endl;
//			std::cout << "Refined Grid size (inner): " << myHestonSolver->getNumberInnerGridPoints() << std::endl;
//		}
//	}
//	else
//	{
//		std::cout << "An unsupported refinement mode has be chosen!" << std::endl;
//		std::cout << "Skipping initial grid refinement!" << std::endl;
//	}

	numGridPoints = myHestonSolver->getNumberGridPoints();

	// Set stochastic data
	myHestonSolver->setStochasticData(theta, kappa, xi, hMatrix, r);

	// Init the grid with on payoff function
	myHestonSolver->initGridWithPayoff(*alpha, dStrike, payoffType);

	// Gridpoints @Money
	std::cout << "Gridpoints @Money: " << myHestonSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;

	// Print interpolation-error at strike
	//	if (numberOfAssets == 2 && coordsType == "cart")
	//	{
	//		myHestonSolver->printPayoffInterpolationError2D(*alpha, "interpolation_error_16384.out", 16384, dStrike);
	//	}

	// Print the payoff function into a gnuplot file
	// Calculate analytic solution if current option is an 1D option
	//	if (numberOfAssets == 1)
	//	{
	//		double minStock = 0.0;
	//		if (payoffType == "std_euro_call")
	//		{
	//			std::vector< std::pair<double, double> > prems;
	//			myHestonSolver->solve1DAnalytic(prems, minStock, maxStock, maxStock/50, dStrike, ((double)(timesteps))*stepsize, true);
	//			myHestonSolver->print1DAnalytic(prems, "analyticHeston.gnuplot");
	//		}
	//		if (payoffType == "std_euro_put")
	//		{
	//			std::vector< std::pair<double, double> > prems;
	//			myHestonSolver->solve1DAnalytic(prems, minStock, maxStock, maxStock/50, dStrike, ((double)(timesteps))*stepsize, false);
	//			myHestonSolver->print1DAnalytic(prems, "analyticHeston.gnuplot");
	//		}
	//	}

	sg::base::DataVector* alphaExact;
	if(numberOfAssets == 1 && payoffType == "std_euro_call")
	{
		//				alphaExact = new sg::base::DataVector(myHestonSolver->getNumberGridPoints());
		//				myHestonSolver->EvaluateHestonExactSurface(*alphaExact,timesteps*stepsize);
		//				myHestonSolver->printGrid(*alphaExact, PLOT_RESOLUTION, "hestonExact.gnuplot");

		//		sg::base::DataVector* alphaCompare = new sg::base::DataVector(myHestonSolver->getNumberGridPoints());
		//		myHestonSolver->CompareHestonBsExact(*alphaCompare, timesteps*stepsize);
		//		myHestonSolver->printGrid(*alphaCompare, 35, "hestonBsCompare.gnuplot");

		//		myHestonSolver->CompareHestonBs1d(timesteps*stepsize, 0.1);
	}

	if (numberOfAssets < 2)
	{
		myHestonSolver->printGrid(*alpha, 100, "payoff.gnuplot");
	}

	if (numberOfAssets < 2)
	{
		myHestonSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
		myHestonSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);

		if (coordsType == "log")
		{
			myHestonSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.gnuplot", true);
			myHestonSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.gnuplot", false);
		}
		if (coordsType == "PAT")
		{
			myHestonSolver->printSparseGridPAT(*alpha, "payoff_surplus_cart.PAT.grid.gnuplot", true);
			myHestonSolver->printSparseGridPAT(*alpha, "payoff_nodal_cart.PAT.grid.gnuplot", false);
		}
	}

	// Start solving the Black Scholes Equation
	if (Solver == "CrNic")
	{
		myHestonSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
	}
	else
	{
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (numberOfAssets < 2)
	{
		// Print the solved Heston Equation into a gnuplot file
		myHestonSolver->printGrid(*alpha, PLOT_RESOLUTION, "solvedHeston.gnuplot");
	}

	//	sg::base::DataVector* alphaCompare = new sg::base::DataVector(myHestonSolver->getNumberGridPoints());
	//	sg::base::DataVector* alphaBsRef = new sg::base::DataVector(myHestonSolver->getNumberGridPoints());
	//	myHestonSolver->CompareHestonNumericToBsExact(*alpha, *alphaBsRef, *alphaCompare, timesteps*stepsize);
	//	myHestonSolver->printGrid(*alphaCompare, PLOT_RESOLUTION, "hestonBsCompare_error.gnuplot");
	//	myHestonSolver->printGrid(*alphaBsRef, PLOT_RESOLUTION, "hestonBsCompare_bsref.gnuplot");


	// Set alphaDone
	//	alphaDone = abs(myHestonSolver->EvalSinglePoint1Asset(sProbe, vProbe, *alpha) - myHestonSolver->EvaluateHestonPriceExact(exp(sProbe), vProbe, xi.get(0), theta.get(0), kappa.get(0), hMatrix.get(0,1), r, timesteps*stepsize, dStrike));
//	alphaDone = myHestonSolver->EvalSinglePoint1Asset(sProbe, vProbe, *alpha);

//	std::cout << "Exact solution: " << myHestonSolver->EvaluateHestonPriceExact(exp(sProbe), vProbe, xi.get(0), theta.get(0), kappa.get(0), hMatrix.get(0,1), r, timesteps*stepsize, dStrike) << std::endl;

	std::stringstream sstm;
	sstm << "solExactDiff" << level << ".gnuplot";
	//	myHestonSolver->CompareHestonSolutionToExact(alpha, alphaExact, sstm.str(), PLOT_RESOLUTION);


	//	if (numberOfAssets == 1 && payoffType == "std_euro_call")
	//	{
	//		// Print the error into a gnuplot file
	//		alpha->sub(*alphaExact);
	//		myHestonSolver->printGrid(*alpha, 50, "hestonError.gnuplot");
	//	}

	if (numberOfAssets < 2)
	{
		myHestonSolver->printSparseGrid(*alpha, "solvedHeston_surplus.grid.gnuplot", true);
		myHestonSolver->printSparseGrid(*alpha, "solvedHeston_nodal.grid.gnuplot", false);

		if (coordsType == "log")
		{
			myHestonSolver->printSparseGridExpTransform(*alpha, "solvedHeston_surplus_cart.grid.gnuplot", true);
			myHestonSolver->printSparseGridExpTransform(*alpha, "solvedHeston_nodal_cart.grid.gnuplot", false);
		}
		if (coordsType == "PAT")
		{
			myHestonSolver->printSparseGridPAT(*alpha, "solvedHeston_surplus_cart.PAT.grid.gnuplot", true);
			myHestonSolver->printSparseGridPAT(*alpha, "solvedHeston_nodal_cart.PAT.grid.gnuplot", false);
		}
	}

	// Test option @ the money
	std::vector<double> point;
	for (size_t i = 0; i < numAssets; i++)
	{
		point.push_back(dStrike);
		double middleVol = (myBoundaries[2*i+1].leftBoundary + myBoundaries[2*i+1].rightBoundary)/2.0;
		point.push_back(middleVol);
	}
	std::cout << "Optionprice at testpoint (Strike): " << myHestonSolver->evalOption(point, *alpha) << std::endl << std::endl;

	//	system("gnuplot /home/sam/Documents/Heston/solExactDiff.cmd");

	delete alpha;
	delete myHestonSolver;
	delete myBoundingBox;
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

			fileStoch.assign(argv[6]);
			fileBound.assign(argv[5]);
			dStrike = atof(argv[7]);
			payoff.assign(argv[8]);
			solver.assign(argv[12]);

			std::string coordsType;
			coordsType.assign(argv[2]);

			// Testing: We'd like to see how the error of a certain point behaves as we increase the region size.
			// For this we keep the stoch file the same and vary the bound file

			// So what I would like to do here is run the BS solver for multiple values of sigma.
			// I'll run it from sigma = 0.3 to sigma = 0.5, which corresponds to the heston v=0.09 to v=0.25
			// I'll collect the results and plot them as a surface for comparison with the Heston results

			vProbe = 0.305;
			sProbe = 0.0;

			const int numTests = 20;
			double initSHalfWidth = 0.6;
			double initVHalfWidth = 0.04;
			double dS = 0.2;
			double dV = 0.05;
			//			double initDiff = 0.005;
			//			double vMins[numTests] = {vProbe - initDiff, vProbe - 2*initDiff, vProbe - 4*initDiff, vProbe - 8*initDiff, vProbe - 16*initDiff, vProbe - 32*initDiff, vProbe - 64*initDiff, vProbe - 128*initDiff};
			//			double vMaxs[numTests] = {vProbe + initDiff, vProbe + 2*initDiff, vProbe + 4*initDiff, vProbe + 8*initDiff, vProbe + 16*initDiff, vProbe + 32*initDiff, vProbe + 64*initDiff, vProbe + 128*initDiff};

			//			initDiff = 0.01;
			//			double sMins[numTests] = {sProbe - initDiff, sProbe - 2*initDiff, sProbe - 4*initDiff, sProbe - 8*initDiff, sProbe - 16*initDiff, sProbe - 32*initDiff, sProbe - 64*initDiff , sProbe - 128*initDiff};
			//			double sMaxs[numTests] = {sProbe + initDiff, sProbe + 2*initDiff, sProbe + 4*initDiff, sProbe + 8*initDiff, sProbe + 16*initDiff, sProbe + 32*initDiff, sProbe + 64*initDiff, sProbe + 128*initDiff};

			//			std::ofstream convFile;
			//			convFile.open("/home/sam/workspace/Heston/convergence.gnuplot");
			//
			//			for(int i=0;i<numTests;i++)
			//			{
			//				std::cout << "Starting test " << i << std::endl;
			//				std::ofstream fileout;
			//				fileout.open("/home/sam/Documents/Heston/tmpBound.bound");
			//				fileout << (sProbe - initSHalfWidth - (i+1)*dS) << " " << (sProbe + initSHalfWidth + (i+1)*dS) << std::endl;
			////				fileout << "-2.04 1.95" << std::endl;
			//				fileout << "0.01 0.61" << std::endl;
			////				fileout << (vProbe - initVHalfWidth - (i+1)*dV) << " " << (vProbe + initVHalfWidth + (i+1)*dV) << std::endl;
			////				fileout << 0.01 << " " << (0.01 + (i+1)*dV) << std::endl;
			//				fileout.close();
			////				vProbe = (0.01 + (0.01 + (i+1)*dV) / 2.0);
			//				testNUnderlyings(atoi(argv[3]), atoi(argv[4]), fileStoch, "/home/sam/Documents/Heston/tmpBound.bound", dStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, coordsType);
			//				convFile << i << " " << alphaDone << std::endl;
			//			}
			//			convFile.close();


			//			std::ofstream convFile;
			//			convFile.open("/home/sam/workspace/Heston/convergence.gnuplot");
			//
			//			for(int i=2;i<10;i++)
			//			{
			//				std::cout << "Starting test " << i << std::endl;
			//				testNUnderlyings(atoi(argv[3]), i, fileStoch, fileBound, dStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, coordsType);
			//				convFile << i << " " << alphaDone << std::endl;
			//			}
			//			convFile.close();

			// Adaptivity tests
//			std::ofstream convFile;
//			convFile.open("/home/sam/workspace/Heston/convergence.gnuplot");
//
//			for(int i=1;i<7;i++)
//			{
//				std::cout << "Starting test " << i << std::endl;
//				refinementThresh = pow(10.0, 0 - i);
//				testNUnderlyings(atoi(argv[3]), atoi(argv[4]), fileStoch, fileBound, dStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, coordsType);
//				convFile << i << " " << numGridPoints << " " << alphaDone << std::endl;
//			}
//			convFile.close();

			testNUnderlyings(atoi(argv[3]), atoi(argv[4]), fileStoch, fileBound, dStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, coordsType);

//			std::cout << "Error: " << alphaDone << std::endl;
		}
	}
	else
	{
		writeHelp();
	}

	return 0;
}
