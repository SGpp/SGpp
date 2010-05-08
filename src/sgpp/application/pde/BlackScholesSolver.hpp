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

#ifndef BLACKSCHOLESSOLVER_HPP
#define BLACKSCHOLESSOLVER_HPP

#include "sgpp.hpp"

#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "grid/type/LinearGrid.hpp"
#include "grid/common/BoundingBox.hpp"

#include "tools/common/StdNormalDistribution.hpp"

#include "application/common/ScreenOutput.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

namespace sg
{

/**
 * This class provides a simple-to-use solver of the multi dimensional Black
 * Scholes Equation that uses Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the Black Scholes
 * Equation with Sparse Grids!
 *
 * @version $HEAD$
 */
class BlackScholesSolver
{
private:
	/// vector that contains the sparse's grid coefficients
	DataVector* alpha;
	/// vector that contains the assets' weight
	DataVector* mus;
	/// vector that contains the standard deviations
	DataVector* sigmas;
	/// Matrix that contains the correlations
	DataVector* rhos;
	/// the riskfree rate
	double r;
	/// the size of one timestep
	double timestepSize;
	/// The number of timesteps that are executed during solving
	size_t nTimesteps;
	/// The Sparse Grid needed in this classificator
	Grid* myGrid;
	/// the number of levels used for an regular grid
	size_t levels;
	/// the dimension of the grid
	size_t dim;
	/// stores if the stochastic asset data was passed to the solver
	bool bStochasticDataAlloc;
	/// stores if the grid was created inside the solver
	bool bGridConstructed;
	/// Stores Pointer to the Grid's Bounding Box
	BoundingBox* myBoundingBox;
	/// Stores Pointer to the Girs's Storage
	GridStorage* myGridStorage;
	/// screen object used in this solver
	ScreenOutput* myScreen;

	/**
	 * returns the option value (payoff value) for an European call option
	 *
	 * @param assetValue the current asset's value
	 * @param strike the strike price of the option
	 *
	 * @return the call premium
	 */
	double get1DEuroCallPayoffValue(double assetValue, double strike);

	/**
	 * This function is a recursive implementation in order the build the evaluation cuboid
	 *
	 * @param evalPoints vector of dynamic size into which the points are "submitted" during calculation
	 * @param curPoint a current point in the d-dimensional space which which is adjusted during this recursive calculations
	 * @param center the center of the cuboid
	 * @param size the precentage of the whole array the cuboid will cover in a every dimension
	 * @param points number of points used in every dimension
	 */
	void getCuboidEvalPoints(std::vector<DataVector>& evalPoints, DataVector& curPoint, std::vector<double>& center, double size, size_t points, size_t curDim);

public:
	/**
	 * Std-Constructor of the solver
	 */
	BlackScholesSolver();

	/**
	 * Std-Destructor of the solver
	 */
	~BlackScholesSolver();

	/**
	 * Use this routine the construct a regular grid to solve the multi-dimensional Black Scholes Equation
	 *
	 * @param myBoundingBox reference to a bounding box that describes the grid
	 * @param level number of the regular's grid levels
	 */
	void constructGrid(BoundingBox& myBoundingBox, size_t level);

	/**
	 * Sets the grid used in this BlackScholes Solver by an given serialized string
	 * of the grid.
	 *
	 * @param serializedGrid a string that describes the grid that should be used in this solver
	 */
	void setGrid(std::string serializedGrid);

	/**
	 * gets the a string the describes the grid which is currently used to solve
	 *
	 * @return a string containing a serialized grid
	 */
	std::string getGrid();

	/**
	 * deletes the grid created within that solver
	 */
	void deleteGrid();

	/**
	 * This function tries to refine the grid such that
	 * most of the grid points are used for interpolation of the singularity. So this grid
	 * is able to approximate the start solution better.
	 *
	 * After refining the grid the payoff function is applied to the grid.
	 *
	 * @param alpha reference to a DataVector object that contains the gird ansatzfunction's coefficients
	 * @param strike pointer to an array the contains all call options strikes
	 * @param payoffType the type of payoff Function used ONLY supported: avgM
	 * @param dStrikeDistance the max. distance from "at the money" a point is allowed to have in order to get refined
	 */
	void refineInitialGrid(DataVector& alpha, double* strike, std::string payoffType, double dStrikeDistance);

	/**
	 * Refines a grid by taking the grid's coefficients into account. This refinement method
	 * refines the grid based on the surplus by refining grid points with big surpluses
	 * first. The number of grid points to refine is specified by a max. percentage
	 * of all grid points.
	 *
	 * @param alpha a DataVector containing the grids coefficients
	 * @param dPercentage percentage of number the give the number of grid points that should be refined
	 */
	void refineInitialGridSurplus(DataVector& alpha, double dPercentage);

	/**
	 * Use this routine the construct a regular grid to solve the multi-dimensional Black Scholes Equation
	 *
	 * Use this routine if you want to solve a problem stored in the format provided by the solving system
	 * released by the University of Bonn, Germany
	 *
	 * @param tfilename absolute path of the file
	 * @param emptyAlpha reference to a DataVector object that contains no elements
	 * @param ishierarchized is set to true if alpha contains surplus after reading the file, otherwise false
	 */
	void constructGrid(std::string tfilename, DataVector& emptyAlpha, bool& ishierarchized);

	/**
	 * Use this routine if you wnat to store a grid in the format provided by the solving system
	 * released by the University of Bonn, Germany
	 *
	 * @param tfilename absolute path of the file
	 * @param alpha reference to a DataVector object that contains the gird ansatzfunction's coefficients
	 * @param ishierarchized set to true, the export is done on the nodal basis
	 */
	void storeGrid(std::string tfilename, DataVector& alpha, bool ishierarchized);

	/**
	 * In order to solve the multi dimensional Black Scholes Equation you have to provided
	 * some statistical data about the underlying (assets' weight, standard deviation
	 * and the correlation between them). This function allows you to set this data.
	 *
	 * @param mus a DataVector that contains the underlyings' weight
	 * @param sigmas a DataVector that contains the underlyings' standard deviations
	 * @param rhos a DataVector that contains the correlations between the underlyings
	 * @param r the riskfree rate used in the market model
	 */
	void setStochasticData(DataVector& mus, DataVector& sigmas, DataVector& rhos, double r);

	/**
	 * Call this routine to use an explicit Euler algorithm to solve the multi dimensional
	 * Black Scholes Equation.
	 *
	 * @param numTimesteps the number of timesteps that should be executed
	 * @param timestepsize the size of the interval one timestep moves forward
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the CG
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param verbose enables verbose output during solving
	 * @param generateAnimation set this to true, if you want to generate a grid output in every timestep
	 * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
	 */
	void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	/**
	 * Call this routine to use an explicit Euler algorithm to solve the multi dimensional
	 * Black Scholes Equation.
	 *
	 * @param numTimesteps the number of timesteps that should be executed
	 * @param timestepsize the size of the interval one timestep moves forward
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the CG
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param verbose enables verbose output during solving
	 * @param generateAnimation set this to true, if you want to generate a grid output in every timestep
	 * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
	 */
	void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	/**
	 * Call this routine to use the Crank Nicolson algorithm to solve the multi dimensional
	 * Black Scholes Equation.
	 *
	 * @param numTimesteps the number of timesteps that should be executed
	 * @param timestepsize the size of the interval one timestep moves forward
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the CG
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 */
	void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha);

	/**
	 * Determines the price of an option in the d-dimensional space
	 *
	 * @param evalPoint coordinates of the point at which the option price should be determined
	 * @param alpha the ansatzfunctions' coefficients
	 *
	 * @return price of option for given point
	 */
	double getOptionPrice(std::vector<double>& evalPoint, DataVector& alpha);

	/**
	 * Evaluates the sparse grid's function given by the stored grid and the alpha coefficients.
	 * on different points specified in EvaluationPoints and stores the result into OptionPrices.
	 *
	 * @param alpha the sparse grid's coefficients
	 * @param OptionPrices DataVector into the which the result of function's evaluation is stored
	 * @param EvaluationPoints DataVector that contains the points at which the sparse grid's function is evaluated
	 */
	void getOptionPricesCuboid(DataVector& alpha, DataVector& OptionPrices, DataVector& EvaluationPoints);

	/**
	 * This function builds an cuboid which will be stored into the EvaluationPoint
	 * variable of this function.
	 * This is by done by building a cuboid around a given center. The size
	 * of the cuboid is determined in every dimension by a fix percent size of the interval in that dimension.
	 * In addition there is a fix number of EvalutionPoints in every dimension. Be aware that this
	 * function returns point to the power of d points.
	 *
	 * @param EvaluationPoints DataVector that will contain the evaluation points afterwards
	 * @param center the center of the cuboid
	 * @param size the precentage of the whole array the cuboid will cover in a every dimension
	 * @param points number of points used in every dimension
	 */
	void getEvaluationCuboid(DataVector& EvaluationPoints, std::vector<double>& center, double size, size_t points);

	/**
	 * Solves the closed form of the Black Scholes equation, the Black Scholes
	 * formular. It evaluates the Black Scholes formular in a Stock Price Range
	 * from 0.0 to maxStock, by increasing the stock price in every step by
	 * a given (small) values, so the analytical solution of the PDE can
	 * be determined and compared.
	 *
	 * @param premiums the result vector, here the combinations of stock price and premium are stored
	 * @param maxStock the maximum stock regarded in these calculations
	 * @param StockInc the increase of the stockprice in one step
	 * @param strike the strike price of the Option
	 * @param t time to maturity
	 */
	void solve1DAnalytic(std::vector< std::pair<double, double> >& premiums, double maxStock, double StockInc, double strike, double t);

	/**
	 * Writes the premiums into a file that can be easily plot with gnuplot
	 *
	 * @param premiums the result vector, here the combinations of stock price and premium are stored
	 * @param tfilename absolute path to file into which the grid's evaluation is written
	 */
	void print1DAnalytic(std::vector< std::pair<double, double> >& premiums, std::string tfilename);

	/**
	 * This is some kind of debug functionality. It writes a file,
	 * that can be used with gnuplot the print the grid.
	 *
	 * Is only implemented for 1D and 2D grids!
	 *
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param PointesPerDimension the distance between evalution points
	 * @param tfilename absolute path to file into which the grid's evaluation is written
	 */
	void printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename);

	/**
	 * Prints the Grid Points of the Sparse Grid either with their node basis value
	 * or their hierarchical surplus
	 *
	 * This function is available for all dimensions
	 *
	 * @param alpha the coefficients of the grid's ansatzfunctions
	 * @param tfilename absoulte path to the file the grid is written into
	 * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
	 */
	void printSparseGrid(DataVector& alpha, std::string tfilename, bool bSurplus);

	/**
	 * Inits the alpha vector with a payoff function of an European call option
	 *
	 * @param alpha the coefficient vector of the grid's ansatzfunctions
	 * @param strike pointer to an array the contains all call options strikes
	 * @param payoffType specifies the type of the combined payoff function; max or avg are available
	 */
	void initGridWithEuroCallPayoff(DataVector& alpha, double* strike, std::string payoffType);

	/**
	 * use this to determine the number of grid points, used to solve
	 * the current problem
	 *
	 * @return the number of grid points
	 */
	size_t getNumberGridPoints();

	/**
	 * use this to determine the number of inner grid points, used to solve
	 * the current problem
	 *
	 * @return the number of inner grid points
	 */
	size_t getNumberInnerGridPoints();

	/**
	 * use this the determine the number of dimensions that are currently used
	 * in the solver.
	 *
	 * @return returns the number of the grid's dimensions, if the grid isn't constructed, yet it returns 0
	 */
	size_t getNumberDimensions();

	/**
	 * Inits the screen object
	 */
	void initScreen();
};

}

#endif /* BLACKSCHOLESSOLVER_HPP */
