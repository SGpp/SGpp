/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HEATEQUATIONSOLVER_HPP
#define HEATEQUATIONSOLVER_HPP

#include "sgpp.hpp"

#include "application/pde/ParabolicPDESolver.hpp"

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
 * Heat Equation on Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the
 * Heat Equation on Sparse Grids!
 *
 * @version $HEAD$
 */
class HeatEquationSolver : public ParabolicPDESolver
{
private:
	/// the heat coefficient
	double a;
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
	virtual ~HeatEquationSolver();

	void constructGrid(BoundingBox& myBoundingBox, size_t level);

	void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul = 0);

	/**
	 * This method sets the heat coefficient of the regarded material
	 *
	 * @param a the heat coefficient
	 */
	void setHeatCoefficient(double a);

	/**
	 * Inits the grid with an smooth heat distribution based the
	 * normal distribution formula
	 *
	 * @param alpha reference to the coefficient's vector
	 * @param mu the exspected value of the normal distribution
	 * @param sigma the sigma of the normal distribution
	 * @param factor a factor that is used to stretch the function values
	 */
	void initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor);

	/**
	 * grid initialization for rotating laser test case
	 *
	 * @param alpha reference to the coefficient's vector
	 * @param heat value of laser's heat
	 * @param nRefinements number of initial refinement before solving the heat equation
	 * @param heat_length expansion of the heat in every dimensions
	 */
	void refineInitialGridWithLaserHeat(DataVector& alpha, double heat, size_t nRefinements, double heat_length);

	/**
	 * Inits the screen object
	 */
	void initScreen();

	/**
	 * Routine to export the matrix of the inner system in matrix
	 * market format
	 *
	 * @param alpha the sparse grid's coefficients
	 * @param tFilename file into which the matrix is written
	 * @param timestepsize the size of the timesteps
	 */
	void storeInnerMatrix(DataVector& alpha, std::string tFilename, double timestepsize);

	/**
	 * Routine to export the matrix of the inner system in matrix
	 * market format
	 *
	 * @param alpha the sparse grid's coefficients
	 * @param tFilename file into which the matrix is written
	 * @param timestepsize the size of the timesteps
	 */
	void storeInnerMatrixDiagonal(DataVector& alpha, std::string tFilename, double timestepsize);

	/**
	 * Routine to export the matrix of the inner system in matrix
	 * market format
	 *
	 * @param alpha the sparse grid's coefficients
	 * @param tFilename file into which the matrix is written
	 * @param timestepsize the size of the timesteps
	 */
	void storeInnerMatrixDiagonalRowSum(DataVector& alpha, std::string tFilename, double timestepsize);

	/**
	 * Routine to export the RHS of the inner system which has to be
	 * solved in order to solve the Poisson equation
	 *
	 * @param alpha the start solution
	 * @param tFilename file into which the rhs is written
	 * @param timestepsize the size of the timesteps
	 */
	void storeInnerRHS(DataVector& alpha, std::string tFilename, double timestepsize);

	/**
	 * Routine to export the solution of the inner system which
	 * has been calculated by Up/Down scheme
	 *
	 * @param alpha the start solution
	 * @param numTimesteps number timesteps
	 * @param timestepsize size of timesteps
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the C
	 * @param tFilename file into which the rhs is written
	 */
	void storeInnerSolution(DataVector& alpha, size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, std::string tFilename);
};

}

#endif /* HEATEQUATIONSOLVER_HPP */
