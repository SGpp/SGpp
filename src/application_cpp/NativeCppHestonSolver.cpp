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

// @todo (heinecke) remove global variables
std::string tFileEvalCuboid = "evalCuboid.data";
std::string tFileEvalCuboidValues = "evalCuboidValues.data";

/// default number of Implicit Euler steps before starting with Crank Nicolson approach
#define CRNIC_IMEUL_STEPS 3
/// default value for epsilon in gridpoints @money
#define DFLT_EPS_AT_MONEY 0.0
/// default value for sigma of refinement normal distribution
#define DFLT_SIGMA_REFINE_NORMDIST 0.15

using namespace sg;
using namespace sg::pde;
using namespace std;

//#define EXPORT_MATRIX_FILES

/**
 * Calls the writeHelp method in the BlackScholesSolver Object
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
//int readStochasticData(std::string tFile, size_t numAssets, sg::base::DataVector& thetas, sg::base::DataVector& volvols, sg::base::DataVector& kappas, sg::base::DataMatrix& rho)
//{
//	std::fstream file;
//	double cur_theta;
//	double cur_kappa;
//	double cur_volvol;
//	double cur_rho;
//
//	file.open(tFile.c_str());
//
//	if(!file.is_open())
//	{
//		std::cout << "Error cannot read file: " << tFile << std::endl;
//		return -1;
//	}
//	// Get number of elements in stoch file, must be numAssests*(numAssests+2)
//	size_t t = 0;
//	double test;
//	do
//	{
//		file >> test;
//		t++;
//	} while (!file.eof());
//	file.close();
//	if (t < ((numAssets*numAssets)+(2*numAssets)))
//	{
//		std::cout << "Invalid stoch file: " << tFile << " Last Value:" << test << std::endl;
//		return -1;
//	}
//
//	file.open(tFile.c_str());
//	for (size_t i = 0; i < numAssets; i++)
//	{
//		file >> cur_theta;
//		file >> cur_kappa;
//		thetas.set(i, cur_theta);
//		kappas.set(i, cur_kappa);
//		for (size_t j = 0; j < numAssets; j++)
//		{
//			file >> cur_rho;
//			rho.set(i,j, cur_rho);
//		}
//	}
//
//	file.close();
//
//	return 0;
//}

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
 * reads the analyze configuration from a file
 *
 * @param tFile the file that contains the analyze data
 * @param numAssets the number of assets
 * @param BoundaryArray Pointer to the Bounding Box array
 * @param points variable to store the number of points in every dimension
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
//int readAnalyzeData(std::string tFile, size_t numAssets, sg::base::DimensionBoundary* BoundaryArray, size_t& points)
//{
//	std::fstream file;
//	double cur_right;
//	double cur_left;
//
//	file.open(tFile.c_str());
//
//	if(!file.is_open())
//	{
//		std::cout << "Error cannot read file: " << tFile << std::endl;
//		return -1;
//	}
//
//	// Get number of elements in analyze file, must be 2*numAssests+1
//	size_t j = 0;
//	double test;
//	do
//	{
//		file >> test;
//		j++;
//	} while (!file.eof());
//	file.close();
//	if (j < (numAssets*2)+1)
//	{
//		std::cout << "Invalid analyze file: " << tFile << " Last Value:" << test << std::endl;
//		return -1;
//	}
//	file.open(tFile.c_str());
//
//
//	for (size_t i = 0; i < numAssets; i++)
//	{
//		file >> cur_left;
//		file >> cur_right;
//
//		BoundaryArray[i].leftBoundary = cur_left;
//		BoundaryArray[i].rightBoundary = cur_right;
//		BoundaryArray[i].bDirichletLeft = true;
//		BoundaryArray[i].bDirichletRight = true;
//	}
//
//	file >> points;
//
//	file.close();
//
//	return 0;
//}

/**
 * reads a cuboid defined by several points from a file. These points are stored in the
 * cuboid DataMatrix
 *
 * @param cuboid DataMatrix into which the evaluations points are stored
 * @param tFile file that contains the cuboid
 * @param dim the dimensions of cuboid
 */
//int readEvalutionCuboid(sg::base::DataMatrix& cuboid, std::string tFile, size_t dim)
//{
//	std::fstream file;
//	double cur_coord;
//
//	file.open(tFile.c_str());
//
//	if(cuboid.getNcols() != dim)
//	{
//		std::cout << "Cuboid-definition file doesn't match: " << tFile << std::endl;
//		return -1;
//	}
//
//	if(!file.is_open())
//	{
//		std::cout << "Error cannot read file: " << tFile << std::endl;
//		return -1;
//	}
//
//	// Get number of lines and resize DataMatrix
//	size_t i = 0;
//	while (!file.eof())
//	{
//		for (size_t d = 0; d < dim; d++)
//		{
//			file >> cur_coord;
//		}
//		i++;
//	}
//	file.close();
//	cuboid.resize(i);
//
//	// Read data from file
//	file.open(tFile.c_str());
//	i = 0;
//	while (!file.eof())
//	{
//		sg::base::DataVector line(dim);
//		line.setAll(0.0);
//		for (size_t d = 0; d < dim; d++)
//		{
//			file >> cur_coord;
//			line.set(d, cur_coord);
//		}
//		cuboid.setRow(i, line);
//		i++;
//	}
//	file.close();
//
//	return 0;
//}

/**
 * reads function values (here option prices) from a file
 *
 * @param values DataVector into which the values will be stored
 * @param tFile file from which the values are read
 * @param numValues number of values stored in the file
 */
//int readOptionsValues(sg::base::DataVector& values, std::string tFile)
//{
//	std::fstream file;
//	double cur_value;
//
//	file.open(tFile.c_str());
//
//	if(!file.is_open())
//	{
//		std::cout << "Error cannot read file: " << tFile << std::endl;
//		return -1;
//	}
//
//	// Count number of lines
//	size_t i = 0;
//	while (!file.eof())
//	{
//		file >> cur_value;
//		i++;
//	}
//	values.resize(i);
//	file.close();
//
//	// Read data from File
//	file.open(tFile.c_str());
//	i = 0;
//	while (!file.eof())
//	{
//		file >> cur_value;
//		values.set(i, cur_value);
//		i++;
//	}
//	file.close();
//
//	return 0;
//}

/**
 * Writes a DataMatrix into a file
 *
 * @param data the DataMatrix that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
//int writeDataMatrix(sg::base::DataMatrix& data, std::string tFile)
//{
//	std::ofstream file;
//	file.open(tFile.c_str());
//
//	if(!file.is_open())
//	{
//		std::cout << "Error cannot write file: " << tFile << std::endl;
//		return -1;
//	}
//
//	for (size_t i = 0; i < data.getNrows(); i++)
//	{
//		for (size_t j = 0; j < data.getNcols(); j++)
//		{
//			file << std::scientific << std::setprecision( 16 ) << data.get(i,j) << " ";
//		}
//		file << std::endl;
//	}
//
//	file.close();
//
//	return 0;
//}


/**
 * Writes a DataVector into a file
 *
 * @param data the DataVector that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
//int writeDataVector(sg::base::DataVector& data, std::string tFile)
//{
//	std::ofstream file;
//	file.open(tFile.c_str());
//
//	if(!file.is_open())
//	{
//		std::cout << "Error cannot write file: " << tFile << std::endl;
//		return -1;
//	}
//
//	for (size_t i = 0; i < data.getSize(); i++)
//	{
//
//		file << std::scientific << std::setprecision( 16 ) << data.get(i) << " " << std::endl;
//	}
//
//	file.close();
//
//	return 0;
//}

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
	//	else if (coordsType == "PAT")
	//	{
	//		myHestonSolver = new sg::finance::HestonSolver(true, true);
	//	}
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

	// init the basis functions' coefficient vector
	sg::base::DataVector* alpha = new sg::base::DataVector(myHestonSolver->getNumberGridPoints());

	std::cout << "Grid has " << level << " Levels" << std::endl;
	std::cout << "Initial Grid size: " << myHestonSolver->getNumberGridPoints() << std::endl;
	std::cout << "Initial Grid size (inner): " << myHestonSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

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
//		alphaExact = new sg::base::DataVector(myHestonSolver->getNumberGridPoints());
//		myHestonSolver->EvaluateHestonExactSurface(*alphaExact,timesteps*stepsize);
//		myHestonSolver->printGrid(*alphaExact, 35, "hestonExact.gnuplot");

//		sg::base::DataVector* alphaCompare = new sg::base::DataVector(myHestonSolver->getNumberGridPoints());
//		myHestonSolver->CompareHestonBsExact(*alphaCompare, timesteps*stepsize);
//		myHestonSolver->printGrid(*alphaCompare, 35, "hestonBsCompare.gnuplot");

//		myHestonSolver->CompareHestonBs1d(timesteps*stepsize, 0.1);
	}

	if (numberOfAssets < 3)
	{
		myHestonSolver->printGrid(*alpha, 100, "payoff.gnuplot");
	}

	if (numberOfAssets < 4)
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
	if (Solver == "ExEul")
	{
		myHestonSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else if (Solver == "ImEul")
	{
		myHestonSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
	}
	else if (Solver == "CrNic")
	{
		myHestonSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
	}
	else if (Solver == "AdBas")
	{
		myHestonSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false);
	}
	else if (Solver == "SCAC")
	{
		myHestonSolver->solveSCAC(timesteps, stepsize, 0.01, CGiterations, CGepsilon, *alpha, false);
	}
	else if (Solver == "SCH")
	{
		myHestonSolver->solveSCH(timesteps, stepsize, 0.0001, CGiterations, CGepsilon, *alpha, false);
	}
	else if (Solver == "SCBDF")
	{
		myHestonSolver->solveSCBDF(timesteps, stepsize, 0.0001, CGiterations, CGepsilon, *alpha, false);
	}
	else if (Solver == "SCEJ")
	{
		myHestonSolver->solveSCEJ(timesteps, stepsize, 0.001, 1.0, CGiterations, CGepsilon, *alpha, false);
	}
	else if (Solver[0] == 'S' && Solver[1] == 'C' && Solver[3] == ':')
	{
		myHestonSolver->solveSC(Solver,timesteps, stepsize, CGiterations, CGepsilon, *alpha, false);
	}
	else
	{
		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
	}

	if (numberOfAssets < 3)
	{
		// Print the solved Heston Equation into a gnuplot file
		myHestonSolver->printGrid(*alpha, 50, "solvedHeston.gnuplot");
	}

//	myHestonSolver->CompareHestonSolutionToExact(alpha, alphaExact, "solExactDiff.gnuplot", 50);


//	if (numberOfAssets == 1 && payoffType == "std_euro_call")
//	{
//		// Print the error into a gnuplot file
//		alpha->sub(*alphaExact);
//		myHestonSolver->printGrid(*alpha, 50, "hestonError.gnuplot");
//	}

	//	if (numberOfAssets < 4)
	//	{
	//		myHestonSolver->printSparseGrid(*alpha, "solvedHeston_surplus.grid.gnuplot", true);
	//		myHestonSolver->printSparseGrid(*alpha, "solvedHeston_nodal.grid.gnuplot", false);
	//
	//		if (coordsType == "log")
	//		{
	//			myHestonSolver->printSparseGridExpTransform(*alpha, "solvedHeston_surplus_cart.grid.gnuplot", true);
	//			myHestonSolver->printSparseGridExpTransform(*alpha, "solvedHeston_nodal_cart.grid.gnuplot", false);
	//		}
	//		if (coordsType == "PAT")
	//		{
	//			myHestonSolver->printSparseGridPAT(*alpha, "solvedHeston_surplus_cart.PAT.grid.gnuplot", true);
	//			myHestonSolver->printSparseGridPAT(*alpha, "solvedHeston_nodal_cart.PAT.grid.gnuplot", false);
	//		}
	//	}

	// Test option @ the money
	std::vector<double> point;
	for (size_t i = 0; i < numAssets; i++)
	{
		point.push_back(dStrike);
		double middleVol = (myBoundaries[2*i+1].leftBoundary + myBoundaries[2*i+1].rightBoundary)/2.0;
		point.push_back(middleVol);
	}
	std::cout << "Optionprice at testpoint (Strike): " << myHestonSolver->evalOption(point, *alpha) << std::endl << std::endl;

	system("gnuplot /home/sam/Documents/Heston/plotHestonSln.cmd");

	delete alpha;
	delete myHestonSolver;
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
 * @param coordsType set the type of coordinates that should be used: cart, log, PAT
 */
//void testNUnderlyingsAnalyze(size_t d, size_t start_l, size_t end_l, std::string fileStoch, std::string fileBound, double dStrike, std::string payoffType,
//		double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, std::string fileAnalyze, std::string coordsType)
//{
//	size_t dim = d;
//	size_t timesteps = timeSt;
//	double stepsize = dt;
//	size_t CGiterations = CGIt;
//	double CGepsilon = CGeps;
//
//	sg::base::DataVector mu(dim);
//	sg::base::DataVector sigma(dim);
//	sg::base::DataMatrix rho(dim,dim);
//
//	sg::base::DataMatrix EvalPoints(1, d);
//
//	double r = riskfree;
//
//	std::vector<sg::base::DataVector> results;
//
//
//	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
//	{
//		return;
//	}
//
//	sg::base::DimensionBoundary* myBoundaries = new sg::base::DimensionBoundary[dim];
//	if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
//	{
//		return;
//	}
//
//	size_t points = 0;
// 	sg::base::DimensionBoundary* myEvalBoundaries = new sg::base::DimensionBoundary[dim];
// 	if (readAnalyzeData(fileAnalyze, dim, myEvalBoundaries, points) != 0)
//	{
//		return;
//	}
//
//	sg::finance::BlackScholesSolver* myBSSolver;
//	if (coordsType == "log")
//	{
//		myBSSolver = new sg::finance::BlackScholesSolver(true);
//	}
//	else if (coordsType == "cart")
//	{
//		myBSSolver = new sg::finance::BlackScholesSolver(false);
//	}
//	else if (coordsType == "PAT")
//	{
//		myBSSolver = new sg::finance::BlackScholesSolver(true, true);
//	}
//	else
//	{
//		// Write Error
//		std::cout << "Unsupported grid transformation!" << std::endl;
//		std::cout << std::endl << std::endl;
//		writeHelp();
//	}
//	sg::base::BoundingBox* myBoundingBox = new sg::base::BoundingBox(dim, myBoundaries);
//	sg::base::BoundingBox* myEvalBoundingBox = new sg::base::BoundingBox(dim, myEvalBoundaries);
//	sg::base::EvalCuboidGenerator* myEvalCuboidGen = new sg::base::EvalCuboidGenerator();
//	delete[] myBoundaries;
//	delete[] myEvalBoundaries;
//
//	// init Screen Object
//	myBSSolver->initScreen();
//
//	for (size_t i = start_l; i <= end_l; i++)
//	{
//		size_t level = i;
//
//		// Reset Solve Time
//		myBSSolver->resetSolveTime();
//
//		// Construct a grid
//		myBSSolver->constructGrid(*myBoundingBox, level);
//
//		// Set stochastic data
//		myBSSolver->setStochasticData(mu, sigma, rho, r);
//
//		// in first iteration -> calculate the evaluation points
//		if (i == start_l)
//		{
//			myEvalCuboidGen->getEvaluationCuboid(EvalPoints, *myEvalBoundingBox, points);
//			writeDataMatrix(EvalPoints, tFileEvalCuboid);
//		}
//
//		// init the basis functions' coefficient vector
//		sg::base::DataVector* alpha = new sg::base::DataVector(myBSSolver->getNumberGridPoints());
//
//		std::cout << "Grid has " << level << " Levels" << std::endl;
//		std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
//		std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;
//
//		// Init the grid with on payoff function
//		myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//
//		// Print the payoff function into a gnuplot file
//		if (dim < 3)
//		{
//			myBSSolver->printGrid(*alpha, 21, "payoff.gnuplot");
//		}
//		if (dim < 4)
//		{
//			myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
//			myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
//
//			if (coordsType == "log")
//			{
//				myBSSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.gnuplot", true);
//				myBSSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.gnuplot", false);
//			}
//			if (coordsType == "PAT")
//			{
//				myBSSolver->printSparseGridPAT(*alpha, "payoff_surplus_cart.PAT.grid.gnuplot", true);
//				myBSSolver->printSparseGridPAT(*alpha, "payoff_nodal_cart.PAT.grid.gnuplot", false);
//			}
//		}
//
//		// Start solving the Black Scholes Equation
//		if (Solver == "ExEul")
//		{
//			myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
//		}
//		else if (Solver == "ImEul")
//		{
//			myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
//		}
//		else if (Solver == "CrNic")
//		{
//			myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
//		}
//		else if (Solver == "AdBas")
//		{
//			myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false);
//		}
//		else if (Solver == "SCAC")
//		{
//			myBSSolver->solveSCAC(timesteps, stepsize, 0.0001, CGiterations, CGepsilon, *alpha, false);
//		}
//		else if (Solver == "SCH")
//		{
//			myBSSolver->solveSCH(timesteps, stepsize, 0.0001, CGiterations, CGepsilon, *alpha, false);
//		}
//		else if (Solver == "SCBDF")
//		{
//			myBSSolver->solveSCBDF(timesteps, stepsize, 0.0001, CGiterations, CGepsilon, *alpha, false);
//		}
//		else if (Solver == "SCEJ")
//		{
//			myBSSolver->solveSCEJ(timesteps, stepsize, 0.001, 1.0, CGiterations, CGepsilon, *alpha, false);
//		}
//		else if (Solver[0] == 'S' && Solver[1] == 'C' && Solver[3] == ':')
//		{
//			myBSSolver->solveSC(Solver,timesteps, stepsize, CGiterations, CGepsilon, *alpha, false);
//		}
//		else
//		{
//			std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
//		}
//
//		if (dim < 3)
//		{
//			// Print the solved Black Scholes Equation into a gnuplot file
//			myBSSolver->printGrid(*alpha, 21, "solvedBS.gnuplot");
//		}
//		if (dim < 4)
//		{
//			myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
//			myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
//
//			if (coordsType == "log")
//			{
//				myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_surplus_cart.grid.gnuplot", true);
//				myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_nodal_cart.grid.gnuplot", false);
//			}
//			if (coordsType == "PAT")
//			{
//				myBSSolver->printSparseGridPAT(*alpha, "solvedBS_surplus_cart.PAT.grid.gnuplot", true);
//				myBSSolver->printSparseGridPAT(*alpha, "solvedBS_nodal_cart.PAT.grid.gnuplot", false);
//			}
//		}
//
//		// Test option @ the money
//		std::vector<double> point;
//		for (size_t x = 0; x < dim; x++)
//		{
//			point.push_back(dStrike);
//		}
//		std::cout << "Optionprice at testpoint (Strike): " << myBSSolver->evalOption(point, *alpha) << std::endl << std::endl;
//
//		if (i == start_l)
//		{
//			// if transformed Black Scholes Equation is used -> transform evaluation domain
//			for (size_t v = 0; v < EvalPoints.getNrows(); v++)
//			{
//				sg::base::DataVector r(dim);
//				EvalPoints.getRow(v, r);
//				myBSSolver->transformPoint(r);
//				EvalPoints.setRow(v, r);
//			}
//		}
//
//		// Evaluate Cuboid
//		sg::base::DataVector Prices(EvalPoints.getNrows());
//		myBSSolver->evaluateCuboid(*alpha, Prices, EvalPoints);
//		results.push_back(Prices);
//
//		// write solution in a additional file
//		std::stringstream level_string;
//		level_string << i;
//		writeDataVector(Prices, tFileEvalCuboidValues+".level_"+ level_string.str());
//		writeDataVector(Prices, tFileEvalCuboidValues);
//
//		if (i > start_l)
//		{
//			std::cout << "=====================================================================" << std::endl;
//			std::cout << "=====================================================================" << std::endl << std::endl;
//			std::cout << "Calculating norms of relative errors to a grid" << std::endl;
//			std::cout << "with " << i << " levels and testing-coboid" << std::endl;
//			std::cout << "with the bounding box:" << std::endl;
//			for (size_t j = 0; j < d; j++)
//			{
//				std::cout << myEvalBoundingBox->getBoundary(j).leftBoundary << " " << myEvalBoundingBox->getBoundary(j).rightBoundary << std::endl;
//			}
//			std::cout << std::endl << std::endl;
//
//			double oldMaxNorm = 0.0;
//			double oldTwoNorm = 0.0;
//
//			std::cout << "Level, max-norm(rel-error), two-norm(rel-error), rate max-norm,  rate two-norm" << std::endl;
//			std::cout << "------------------------------------------------------------------------------" << std::endl;
//			// Calculate relative errors and some norms
//			for (size_t j = 0; j < i-start_l; j++)
//			{
//				sg::base::DataVector maxLevel(results[i-start_l]);
//				sg::base::DataVector relError(results[j]);
//				double maxNorm = 0.0;
//				double l2Norm = 0.0;
//
//				// calculate relative error
//				relError.sub(maxLevel);
//				relError.componentwise_div(maxLevel);
//
//				// calculate max. norm of relative error
//				maxNorm = relError.maxNorm();
//
//				// calculate two norm of relative error
//				l2Norm = relError.RMSNorm();
//
//				// Printing norms
//				//std::cout << "Level " << j + start_l << ": max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << l2Norm << "; rate max-norm: " << log(oldMaxNorm/maxNorm) << "; rate two-norm: " << log(oldTwoNorm/l2Norm) << std::endl;
//				std::cout << j + start_l << ", " << maxNorm << ", " << l2Norm << ", " << log(oldMaxNorm/maxNorm) << ", " << log(oldTwoNorm/l2Norm) << std::endl;
//
//				oldMaxNorm = maxNorm;
//				oldTwoNorm = l2Norm;
//			}
//		}
//		std::cout << std::endl << std::endl;
//
//		myBSSolver->deleteGrid();
//		delete alpha;
//
//		std::cout << std::endl;
//	}
//
//	delete myEvalBoundingBox;
//	delete myEvalCuboidGen;
//	delete myBSSolver;
//	delete myBoundingBox;
//}


/**
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call option, with Initial
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
 * @param refinementMode the mode selected for surplus refinement: available: classic, maxLevel
 * @param maxRefineLevel ignored for refinement mode classic, in maxLevel: max. level to which the grid is refined
 * @param numRefinePoints number of points that should be refined in each refine iteration before Black Scholes Equation is solved: -1 try to refine all points steered by threshold
 * @param nIterAdaptSteps number of the iterative Grid Refinement that should be executed
 * @param dRefineThreshold Threshold for a point's surplus for refining this point
 * @param useCoarsen specifies if the grid should be coarsened between timesteps
 * @param adaptSolvingMode specifies which adaptive methods are applied during solving the BS Equation
 * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
 * @param coordsType set the type of coordinates that should be used: cart, log, PAT
 * @param useNormalDist enable local initial refinement based on a normal distribution
 */
//void testNUnderlyingsAdaptSurplus(size_t d, size_t l, std::string fileStoch, std::string fileBound, double dStrike,
//		std::string payoffType, double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps,
//		std::string Solver, std::string refinementMode, int numRefinePoints, size_t maxRefineLevel, size_t nIterAdaptSteps, double dRefineThreshold,
//		bool useCoarsen, std::string adaptSolvingMode, double coarsenThreshold, std::string coordsType, bool useNormalDist)
//{
//	size_t dim = d;
//	size_t level = l;
//	size_t timesteps = timeSt;
//	double stepsize = dt;
//	size_t CGiterations = CGIt;
//	double CGepsilon = CGeps;
//
//	sg::base::DataVector mu(dim);
//	sg::base::DataVector sigma(dim);
//	sg::base::DataMatrix rho(dim,dim);
//
//	double r = riskfree;
//
//	if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0)
//	{
//		return;
//	}
//
//	sg::base::DimensionBoundary* myBoundaries = new sg::base::DimensionBoundary[dim];
//	if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0)
//	{
//		return;
//	}
//
//	sg::finance::BlackScholesSolver* myBSSolver;
//	if (coordsType == "log")
//	{
//		myBSSolver = new sg::finance::BlackScholesSolver(true);
//	}
//	else if (coordsType == "cart")
//	{
//		myBSSolver = new sg::finance::BlackScholesSolver(false);
//	}
//	else if (coordsType == "PAT")
//	{
//		myBSSolver = new sg::finance::BlackScholesSolver(true, true);
//	}
//	else
//	{
//		// Write Error
//		std::cout << "Unsupported grid transformation!" << std::endl;
//		std::cout << std::endl << std::endl;
//		writeHelp();
//	}
//
//	sg::base::BoundingBox* myBoundingBox = new sg::base::BoundingBox(dim, myBoundaries);
//	delete[] myBoundaries;
//
//	// init Screen Object
//	myBSSolver->initScreen();
//
//	// Construct a grid
//	myBSSolver->constructGrid(*myBoundingBox, level);
//
//	// Enable Coarsening
//	if (useCoarsen == true)
//	{
//		myBSSolver->setEnableCoarseningData(adaptSolvingMode, refinementMode, maxRefineLevel, -1, coarsenThreshold, dRefineThreshold);
//	}
//
//	// init the basis functions' coefficient vector
//	sg::base::DataVector* alpha = new sg::base::DataVector(myBSSolver->getNumberGridPoints());
//
//	// Set stochastic data
//	myBSSolver->setStochasticData(mu, sigma, rho, r);
//
//	// Init the grid with on payoff function
//	myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//
//	std::vector<double> norm_mu;
//	std::vector<double> norm_sigma;
//	double refineSigma = DFLT_SIGMA_REFINE_NORMDIST;
//
//	// estimate refine sigma from evaluation cuboid
//	// read Evaluation cuboid
//	sg::base::DataMatrix EvalCuboid(1, dim);
//	int retCuboid = readEvalutionCuboid(EvalCuboid, tFileEvalCuboid, dim);
//
//	// read reference values for evaluation cuboid
//	sg::base::DataVector EvalCuboidValues(1);
//	int retCuboidValues = readOptionsValues(EvalCuboidValues, tFileEvalCuboidValues);
//
//	if (EvalCuboid.getNrows() != EvalCuboidValues.getSize())
//	{
//		retCuboid = 1;
//		retCuboidValues = 1;
//	}
//
//	if (retCuboid == 0 && retCuboidValues == 0)
//	{
//		refineSigma = EvalCuboid.max(0) - EvalCuboid.min(0);
//	}
//
//	if (useNormalDist == true)
//	{
//		for (size_t i = 0; i < d; i++)
//		{
//			norm_mu.push_back(dStrike);
//			norm_sigma.push_back(refineSigma);
//		}
//	}
//
////	if (useCoarsen == true)
////	{
////		for (size_t i = 0 ; i < nIterAdaptSteps; i++)
////		{
////			myBSSolver->coarsenInitialGridSurplus(*alpha, coarsenThreshold);
////		}
////	}
//
//	std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
//	std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;
//
//	// refine the grid to approximate the singularity in the start solution better
//	if (refinementMode == "classic")
//	{
//		for (size_t i = 0 ; i < nIterAdaptSteps; i++)
//		{
//			std::cout << "Refining Grid..." << std::endl;
//			if (useNormalDist == true)
//			{
//				myBSSolver->refineInitialGridSurplusSubDomain(*alpha, numRefinePoints, dRefineThreshold, norm_mu, norm_sigma);
//			}
//			else
//			{
//				myBSSolver->refineInitialGridSurplus(*alpha, numRefinePoints, dRefineThreshold);
//			}
//			myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//			std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
//			std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
//		}
//
//	}
//	else if (refinementMode == "maxLevel")
//	{
//		size_t oldGridSize = 0;
//		size_t newGridSize = myBSSolver->getNumberGridPoints();
//		size_t addedGridPoint = 0;
//		size_t stepCounter = 0;
//		if (nIterAdaptSteps > 0)
//		{
//			do
//			{
//				oldGridSize = newGridSize;
//				std::cout << "Refining Grid..." << std::endl;
//				if (useNormalDist == true)
//				{
//					myBSSolver->refineInitialGridSurplusToMaxLevelSubDomain(*alpha, dRefineThreshold, maxRefineLevel, norm_mu, norm_sigma);
//				}
//				else
//				{
//					myBSSolver->refineInitialGridSurplusToMaxLevel(*alpha, dRefineThreshold, maxRefineLevel);
//				}
//				myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//				std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
//				std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
//				newGridSize = myBSSolver->getNumberGridPoints();
//				addedGridPoint = newGridSize - oldGridSize;
//				stepCounter++;
//			} while ((addedGridPoint > 0) && (stepCounter < nIterAdaptSteps));
//		}
//	}
//	else
//	{
//		std::cout << "An unsupported refinement mode has be chosen!" << std::endl;
//		std::cout << "Skipping initial grid refinement!" << std::endl;
//	}
//	std::cout << std::endl << std::endl << std::endl;
//
////	if (useCoarsen == true)
////	{
////		myBSSolver->coarsenInitialGridSurplus(*alpha, coarsenThreshold);
////		std::cout << "Coarsened Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
////		std::cout << "Coarsened Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
////		std::cout << std::endl << std::endl << std::endl;
////	}
//
//	// Print the payoff function into a gnuplot file
//	if (dim < 3)
//	{
//		myBSSolver->printGrid(*alpha, 21, "payoff.gnuplot");
//	}
//	if (dim < 4)
//	{
//		myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
//		myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);
//
//		if (coordsType == "log")
//		{
//			myBSSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.gnuplot", true);
//			myBSSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.gnuplot", false);
//		}
//		if (coordsType == "PAT")
//		{
//			myBSSolver->printSparseGridPAT(*alpha, "payoff_surplus_cart.PAT.grid.gnuplot", true);
//			myBSSolver->printSparseGridPAT(*alpha, "payoff_nodal_cart.PAT.grid.gnuplot", false);
//		}
//	}
//
//	// Gridpoints @Money
//	std::cout << "Gridpoints @Money: " << myBSSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;
//
//	// Print interpolation-error at strike
//	if (dim == 2 && coordsType == "log")
//	{
//		myBSSolver->printPayoffInterpolationError2D(*alpha, "interpolation_error_16384.out", 16384, dStrike);
//	}
//
//	// Start solving the Black Scholes Equation
//	if (Solver == "ExEul")
//	{
//		myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
//	}
//	else if (Solver == "ImEul")
//	{
//		myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
//	}
//	else if (Solver == "CrNic")
//	{
//		myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
//	}
//	else if (Solver == "AdBas")
//	{
//		myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false);
//	}
//	else if (Solver == "SCAC")
//	{
//		myBSSolver->solveSCAC(timesteps, stepsize, 0.0001, CGiterations, CGepsilon, *alpha, false);
//	}
//	else if (Solver == "SCH")
//	{
//		myBSSolver->solveSCH(timesteps, stepsize, 0.0001, CGiterations, CGepsilon, *alpha, false);
//	}
//	else if (Solver == "SCBDF")
//	{
//		myBSSolver->solveSCBDF(timesteps, stepsize, 0.0001, CGiterations, CGepsilon, *alpha, false);
//	}
//	else if (Solver == "SCEJ")
//	{
//		myBSSolver->solveSCEJ(timesteps, stepsize, 0.001, 1.0, CGiterations, CGepsilon, *alpha, false);
//	}
//	else if (Solver[0] == 'S' && Solver[1] == 'C' && Solver[3] == ':')
//	{
//		myBSSolver->solveSC(Solver,timesteps, stepsize, CGiterations, CGepsilon, *alpha, false);
//	}
//	else
//	{
//		std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
//	}
//
//	if (dim < 3)
//	{
//		// Print the solved Black Scholes Equation into a gnuplot file
//		myBSSolver->printGrid(*alpha, 21, "solvedBS.gnuplot");
//	}
//	if (dim < 4)
//	{
//		myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.gnuplot", true);
//		myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.gnuplot", false);
//
//		if (coordsType == "log")
//		{
//			myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_surplus_cart.grid.gnuplot", true);
//			myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_nodal_cart.grid.gnuplot", false);
//		}
//		if (coordsType == "PAT")
//		{
//			myBSSolver->printSparseGridPAT(*alpha, "solvedBS_surplus_cart.PAT.grid.gnuplot", true);
//			myBSSolver->printSparseGridPAT(*alpha, "solvedBS_nodal_cart.PAT.grid.gnuplot", false);
//		}
//	}
//
//	// Test option @ the money
//	std::vector<double> point;
//	for (size_t i = 0; i < d; i++)
//	{
//		point.push_back(dStrike);
//	}
//	std::cout << "Optionprice at testpoint (Strike): " << myBSSolver->evalOption(point, *alpha) << std::endl << std::endl;
//
//	// calculate relative errors
//	////////////////////////////
//	double maxNorm = 0.0;
//	double l2Norm = 0.0;
//
//	if (retCuboid == 0 && retCuboidValues == 0)
//	{
//		// If the transformed Black Scholes Equation is used -> transform Eval-cuboid
//		for (size_t v = 0; v < EvalCuboid.getNrows(); v++)
//		{
//			sg::base::DataVector r(dim);
//			EvalCuboid.getRow(v, r);
//			myBSSolver->transformPoint(r);
//			EvalCuboid.setRow(v, r);
//		}
//
//		std::cout << "Calculating relative errors..." << std::endl;
//		// Evaluate Cuboid
//		sg::base::DataVector Prices(EvalCuboid.getNrows());
//		myBSSolver->evaluateCuboid(*alpha, Prices, EvalCuboid);
//
//		sg::base::DataVector relError(Prices);
//
//		// calculate relative error
//		relError.sub(EvalCuboidValues);
//		relError.componentwise_div(EvalCuboidValues);
//
//		// calculate max. norm of relative error
//		maxNorm = relError.maxNorm();
//
//		// calculate two norm of relative error
//		l2Norm = relError.RMSNorm();
//
//		// Printing norms
//		std::cout << "Results: max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << l2Norm << std::endl;
//
//		// reprint data with prefix -> can be easily grep-ed
//		std::cout << std::endl << std::endl;
//	}
//	else
//	{
//		std::cout << "Couldn't open evaluation cuboid data -> skipping tests!" << std::endl << std::endl;
//	}
//
//	std::cout << "$ Startlevel: " << level << "; RefineMode: " << refinementMode << "; MaxRefLevel: " << maxRefineLevel << std::endl;
//	std::string normDistrefine;
//	if (useNormalDist == true)
//	{
//		std::stringstream normDistRefineStream;
//		normDistRefineStream << "solveNDadaptSurplusSubDomain;" << dStrike << ";" << refineSigma;
//		normDistrefine = normDistRefineStream.str();
//		std::cout << "$ AdaptSurplus-Mode: solveNDadaptSurplusSubDomain" << std::endl;
//		std::cout << "$ Refine mu = " << dStrike << "; Refine sigma = " << refineSigma << std::endl;
//	}
//	else
//	{
//		normDistrefine = "solveNDadaptSurplus;-1.0;1.0";
//		std::cout << "$ AdaptSurplus-Mode: solveNDadaptSurplus" << std::endl;
//	}
//
//	std::cout << "$ NumRefinements: " << nIterAdaptSteps << "; RefineThreshd: " << dRefineThreshold << std::endl;
//	std::cout << "$ AdpatSolveMode: " << adaptSolvingMode << "; CoarsenThreshd: " << coarsenThreshold << std::endl;
//	std::cout << "$ Start #gridpoints (inner): " << myBSSolver->getStartInnerGridSize() << std::endl;
//	std::cout << "$ Final #gridpoints (inner): " << myBSSolver->getFinalInnerGridSize() << std::endl;
//	std::cout << "$ Average #gridpoints (inner): " << myBSSolver->getAverageInnerGridSize() << std::endl;
//	std::cout << "$ Needed iterations: " << myBSSolver->getNeededIterationsToSolve() << "; Needed time: " << myBSSolver->getNeededTimeToSolve() << std::endl;
//	std::cout << "$ Results: max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << l2Norm << std::endl;
//	std::cout << "$ Optionprice at testpoint (Strike): " << myBSSolver->evalOption(point, *alpha) << std::endl;
//	std::cout << "$ CSV-DATA: " << level << ";" << refinementMode << ";" << maxRefineLevel << ";" << nIterAdaptSteps
//		<< ";" << dRefineThreshold << ";" << normDistrefine << ";" << adaptSolvingMode << ";" << coarsenThreshold
//		<< ";" << myBSSolver->getStartInnerGridSize() << ";" << myBSSolver->getFinalInnerGridSize()
//		<< ";" << myBSSolver->getAverageInnerGridSize() << ";" << myBSSolver->getNeededIterationsToSolve()
//		<< ";" << myBSSolver->getNeededTimeToSolve() << ";" << maxNorm << ";" << l2Norm << std::endl;
//	std::cout << std::endl << std::endl;
//
//#ifdef EXPORT_MATRIX_FILES
//	// print inner matrix
//	std::stringstream mtxFile;
//	mtxFile.clear();
//	alpha->setAll(0.0);
//	myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//
//	mtxFile << "SG_BlackScholes_InnerMatrix_" << dim << "d_" << level << "l.mtx";
//	myBSSolver->storeInnerMatrix(*alpha, mtxFile.str(), dt);
//
//	// print inner matrix, diagonal
//	std::stringstream mtxFileDiagonal;
//	mtxFileDiagonal.clear();
//	alpha->setAll(0.0);
//	myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//
//	mtxFileDiagonal << "SG_BlackScholes_InnerMatrixDiagonal_" << dim << "d_" << level << "l.mtx";
//	myBSSolver->storeInnerMatrixDiagonal(*alpha, mtxFileDiagonal.str(), dt);
//
//	// print inner matrix, diagonal row sum
//	std::stringstream mtxFileDiagonalRowSum;
//	mtxFileDiagonalRowSum.clear();
//	alpha->setAll(0.0);
//	myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//
//	mtxFileDiagonalRowSum << "SG_BlackScholes_InnerMatrixDiagonalRowSum_" << dim << "d_" << level << "l.mtx";
//	myBSSolver->storeInnerMatrixDiagonalRowSum(*alpha, mtxFileDiagonalRowSum.str(), dt);
//
//	// print inner rhs
//	std::stringstream rhsFile;
//	rhsFile.clear();
//	alpha->setAll(0.0);
//	myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//
//	rhsFile << "SG_BlackScholes_InnerRHS_" << dim << "d_" << level << "l.vec";
//	myBSSolver->storeInnerRHS(*alpha, rhsFile.str(), dt);
//
//	// print inner solution
//	std::stringstream solFile;
//	solFile.clear();
//	alpha->setAll(0.0);
//	myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
//
//	solFile << "SG_BlackScholes_InnerSolution_" << dim << "d_" << level << "l.vec";
//	myBSSolver->storeInnerSolution(*alpha, timesteps, stepsize, CGiterations, CGepsilon, solFile.str());
//	std::cout << std::endl << std::endl;
//#endif
//
//	delete myBSSolver;
//	delete myBoundingBox;
//	delete alpha;
//}

// function to interpolate
double f(int dim, double *x, void *clientdata) {
	double res = 1.0;
	for (int i=0; i<dim; i++) {
		res *= 4*x[i]*(1-x[i]);
	}
	return res;
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

			testNUnderlyings(atoi(argv[3]), atoi(argv[4]), fileStoch, fileBound, dStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, coordsType);
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

			fileStoch.assign(argv[7]);
			fileBound.assign(argv[6]);
			dStrike = atof(argv[8]);
			fileAnalyze.assign(argv[16]);
			payoff.assign(argv[9]);
			solver.assign(argv[13]);

			std::string coordsType;
			coordsType.assign(argv[2]);

			//			testNUnderlyingsAnalyze(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), fileStoch, fileBound, dStrike, payoff, atof(argv[10]), (size_t)(atof(argv[11])/atof(argv[12])), atof(argv[12]), atoi(argv[14]), atof(argv[15]), solver, fileAnalyze, coordsType);
		}
	}
	else if (option == "solveNDadaptSurplus" || option == "solveNDadaptSurplusSubDomain")
	{
		if (argc != 21)
		{
			writeHelp();
		}
		else
		{
			bool isNormalDist = false;
			if (option == "solveNDadaptSurplusSubDomain")
			{
				isNormalDist = true;
			}
			std::string fileStoch;
			std::string fileBound;
			double dStrike;
			std::string ani;
			std::string solver;
			std::string payoff;
			std::string refinementMode;
			std::string adaptSolveMode;

			fileStoch.assign(argv[6]);
			fileBound.assign(argv[5]);
			dStrike = atof(argv[7]);
			payoff.assign(argv[8]);
			solver.assign(argv[12]);
			refinementMode.assign(argv[15]);
			adaptSolveMode.assign(argv[19]);

			std::string coordsType;
			coordsType.assign(argv[2]);

			if (refinementMode != "maxLevel" && refinementMode != "classic")
			{
				std::cout << "Unsupported refinement type! classic or maxLevel are supported!" << std::endl;
				std::cout << std::endl << std::endl;
				writeHelp();
				return 0;
			}

			bool useAdaptSolve = false;
			if (adaptSolveMode == "coarsen" || adaptSolveMode == "refine" || adaptSolveMode == "coarsenNrefine")
			{
				useAdaptSolve = true;
			}
			else if (adaptSolveMode == "none")
			{
				useAdaptSolve = false;
			}
			else
			{
				std::cout << "Unsupported adapt solve mode! none, coarsen, refine or coarsenNrefine are supported!" << std::endl;
				std::cout << std::endl << std::endl;
				writeHelp();
				return 0;
			}

			//			testNUnderlyingsAdaptSurplus(atoi(argv[3]), atoi(argv[4]), fileStoch, fileBound, dStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10])/atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, refinementMode, -1, atoi(argv[16]), atoi(argv[17]), atof(argv[18]), useAdaptSolve, adaptSolveMode, atof(argv[20]), coordsType, isNormalDist);
		}
	}
	else
	{
		writeHelp();
	}

	return 0;
}
