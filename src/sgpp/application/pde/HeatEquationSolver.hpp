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

#ifndef HEATEQUATIONSOLVER_HPP
#define HEATEQUATIONSOLVER_HPP

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
 * This class provides a simple-to-use solver of the multi dimensional
 * Heat Equation that uses Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the
 * Heat Equation with Sparse Grids!
 *
 * @version $HEAD$
 */
class HeatEquationSolver
{
private:
	/// vector that contains the sparse's grid coefficients
	DataVector* alpha;
	/// the heat coefficient
	double a;
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
	/// stores if the grid was created inside the solver
	bool bGridConstructed;
	/// Stores Pointer to the Grid's Bounding Box
	BoundingBox* myBoundingBox;
	/// Stores Pointer to the Girs's Storage
	GridStorage* myGridStorage;
	/// screen object used in this solver
	ScreenOutput* myScreen;

public:
	/**
	 * Std-Constructor of the solver
	 */
	HeatEquationSolver();

	/**
	 * Std-Destructor of the solver
	 */
	~HeatEquationSolver();

	/**
	 * Use this routine the construct a regular grid to solve the multi-dimensional Black Scholes Equation
	 *
	 * @param myBoundingBox reference to a bounding box that describes the grid
	 * @param level number of the regular's grid levels
	 * @param useBoundary use a grid with explicit boundary conditions
	 */
	void constructGrid(BoundingBox& myBoundingBox, size_t level, bool useBoundary);

	/**
	 * Call this routine to use an explicit Euler algorithm to solve the multi dimensional
	 * Heat Equation.
	 *
	 * @param numTimesteps the number of timesteps that should be executed
	 * @param timestepsize the size of the interval one timestep moves forward
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the CG
	 * @param a the heat coefficient of the analyzed material
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param verbose enables verbose output during solving
	 * @param generateAnimation set this to true, if you want to generate a grid output in every timestep
	 * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
	 */
	void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	/**
	 * Call this routine to use an implicit Euler algorithm to solve the multi dimensional
	 * Heat Equation.
	 *
	 * @param numTimesteps the number of timesteps that should be executed
	 * @param timestepsize the size of the interval one timestep moves forward
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the CG
	 * @param a the heat coefficient of the analyzed material
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param verbose enables verbose output during solving
	 * @param generateAnimation set this to true, if you want to generate a grid output in every timestep
	 * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
	 */
	void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	/**
	 * Call this routine to use the Crank Nicolson algorithm to solve the multi dimensional
	 * Heat Equation.
	 *
	 * @param numTimesteps the number of timesteps that should be executed
	 * @param timestepsize the size of the interval one timestep moves forward
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the CG
	 * @param a the heat coefficient of the analyzed material
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 */
	void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha);

	/**
	 * This is some kind of debug functionality. It writes a file,
	 * that can be used with gnuplot the print the grid.
	 *
	 * Is only implemented for 1D and 2D grids!
	 *
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param PointesPerDimension number of points that should be evaluated in one dimension
	 * @param tfilename absolute path to file into which the grid's evaluation is written
	 */
	void printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename);

	/**
	 * use this to determine the number of grid points, used to solve
	 * the current problem
	 *
	 * @return the number of grid points
	 */
	size_t getNumberGridPoints();

	/**
	 * Inits the grid in the middle of the whole domain with one single heat
	 *
	 * @param alpha reference to the coefficients vector
	 * @param heat the value of the heat in the middle of the domain
	 */
	void initGridWithSingleHeat(DataVector& alpha, double heat);

	/**
	 * Inits the grid in the middle the domain with an smooth heat distribution that the
	 * normal distribution formula
	 *
	 * @param alpha reference to the coefficients vector
	 * @param mu the exspected value of the normal distribution
	 * @param sigma the sigma of the normal distribution
	 * @param factor a factor that is used to stretch the function values
	 */
	void initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor);

	/**
	 * Inits the grid with a constant heat
	 *
	 * @param alpha reference to the coefficients vector
	 * @param constHeat the temperature of the constant heat
	 */
	void initGridWithConstantHeat(DataVector& alpha, double constHeat);

	/**
	 * Inits the screen object
	 */
	void initScreen();
};

}

#endif /* HEATEQUATIONSOLVER_HPP */
