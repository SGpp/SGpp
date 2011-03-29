/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LASERHEATEQUATIONSOLVER_HPP
#define LASERHEATEQUATIONSOLVER_HPP

#include "application/pde/HeatEquationSolver.hpp"

namespace sg
{

/**
 * This class provides a simple-to-use solver of the multi dimensional
 * Heat Equation on Sparse Grids.
 *
 * This class provides an extension to the classic HeatEquationSolver.
 * Here the testcase of a rotating laser beam is solved using the
 * heat equation.
 *
 * @version $HEAD$
 */
class LaserHeatEquationSolver : public HeatEquationSolver
{
public:
	/**
	 * Std-Constructor of the solver
	 */
	LaserHeatEquationSolver();

	/**
	 * Std-Destructor of the solver
	 */
	virtual ~LaserHeatEquationSolver();

	void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul = 0);

	/**
	 * grid initialization for rotating laser test case
	 *
	 * @param alpha reference to the coefficient's vector
	 * @param heat value of laser's heat
	 * @param nRefinements number of initial refinement before solving the heat equation
	 * @param maxLevel max. Level in the refined grid
	 * @param heat_sigma expansion of the heat in every dimensions (sigma of normal distribution)
	 */
	void refineInitialGridWithLaserHeat(DataVector& alpha, double heat, size_t nRefinements, size_t maxLevel, double heat_sigma);

	void initScreen();
};

}

#endif /* LASERHEATEQUATIONSOLVER_HPP */
