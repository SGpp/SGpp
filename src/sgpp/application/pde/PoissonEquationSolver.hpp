/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef POISSONEQUATIONSOLVER_HPP
#define POISSONEQUATIONSOLVER_HPP

#include "sgpp.hpp"

#include "application/pde/EllipticPDESolver.hpp"

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
 * Poisson Equation on Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the
 * Poisson Equation on Sparse Grids!
 *
 * @version $HEAD$
 */
class PoissonEquationSolver : public EllipticPDESolver
{
private:
	/// screen object used in this solver
	ScreenOutput* myScreen;

public:
	/**
	 * Std-Constructor of the solver
	 */
	PoissonEquationSolver();

	/**
	 * Std-Destructor of the solver
	 */
	virtual ~PoissonEquationSolver();

	void constructGrid(BoundingBox& myBoundingBox, size_t level);

	void solvePDE(DataVector& alpha, DataVector& rhs, size_t maxCGIterations, double epsilonCG, bool verbose = false);

	/**
	 * Inits the grid with a smooth heat distribution (based on
	 * a std-normal distribution) on its boundaries
	 *
	 * Coefficients of inner grid points are set to zero
	 * since an elliptical PDE is solved
	 *
	 * @param alpha reference to the coefficients vector
	 * @param mu the exspected value of the normal distribution
	 * @param sigma the sigma of the normal distribution
	 * @param factor a factor that is used to stretch the function values
	 */
	void initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor);

	/**
	 * Inits the grid with a smooth heat distribution (based on
	 * a std-normal distribution) on its boundaries
	 *
	 * Coefficients of inner grid points aren't set to zero
	 * since they are used to hint an adaptive refinement
	 * of the grid BEFORE solving the PDE.
	 *
	 * @param alpha reference to the coefficients vector
	 * @param mu the exspected value of the normal distribution
	 * @param sigma the sigma of the normal distribution
	 * @param factor a factor that is used to stretch the function values
	 */
	void initGridWithSmoothHeatFullDomain(DataVector& alpha, double mu, double sigma, double factor);

	/**
	 * Inits the grid with a heat distribution based on
	 * the e-function
	 *
	 * The e-function is shifted in that way the right boundary
	 * values becomes 1 (in case of factor = 1)
	 *
	 * @param alpha reference to the coefficient's vector
	 * @param factor a constant factor used to enlarge the exp-functions input parameter
	 */
	void initGridWithExpHeat(DataVector& alpha, double factor = 1.0);

	/**
	 * Inits the grid with a heat distribution based on
	 * the e-function
	 *
	 * The e-function is shifted in that way the right boundary
	 * values becomes 1 (in case of factor = 1)
	 *
	 * @param alpha reference to the coefficient's vector
	 * @param factor a constant factor used to enlarge the exp-functions input parameter
	 */
	void initGridWithExpHeatFullDomain(DataVector& alpha, double factor = 1.0);

	/**
	 * Routine to export the matrix of the inner system in matrix
	 * market format
	 *
	 * @param tFilename file into which the matrix is written
	 */
	void storeInnerMatrix(std::string tFilename);

	/**
	 * Routine to export inner system matrix's diagonal as
	 * diagonal matrix in matrix market format
	 *
	 * @param tFilename file into which the matrix is written
	 */
	void storeInnerMatrixDiagonal(std::string tFilename);

	/**
	 * Routine to export inner system matrix's row sum as
	 * diagonal matrix in matrix market format
	 *
	 * @param tFilename file into which the matrix is written
	 */
	void storeInnerMatrixDiagonalRowSum(std::string tFilename);

	/**
	 * Routine to export the RHS of the inner system which has to be
	 * solved in order to solve the Poisson equation
	 *
	 * @param alpha the start solution
	 * @param tFilename file into which the rhs is written
	 */
	void storeInnerRHS(DataVector& alpha, std::string tFilename);

	/**
	 * Routine to export the solution of the inner system which
	 * has been calculated by Up/Down scheme
	 *
	 * @param alpha the start solution
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the C
	 * @param tFilename file into which the rhs is written
	 */
	void storeInnerSolution(DataVector& alpha, size_t maxCGIterations, double epsilonCG, std::string tFilename);

	/**
	 * Inits the screen object
	 */
	void initScreen();
};

}

#endif /* POISSONEQUATIONSOLVER_HPP */
