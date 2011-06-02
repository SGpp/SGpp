/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LASERHEATEQUATIONSOLVER2D_HPP
#define LASERHEATEQUATIONSOLVER2D_HPP

#include "application/pde/HeatEquationSolver.hpp"

namespace sg
{
namespace pde
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
class LaserHeatEquationSolver2D : public HeatEquationSolver
{
private:
	/// stores the beam velocity
	double beam_velocity_;
	/// stores the expansions of the laser's heat
	double heat_sigma_;
	/// stores the grid's max. refinement level
	size_t max_level_;
	/// threshold for refineing the grid during solution process
	double refine_threshold_;
	/// threshold for coarsening the grid during solution process
	double coarsen_threshold_;
	/// heating factor
	double heat_;

public:
	/**
	 * Std-Constructor of the solver
	 *
	 * @param beam_velocity the velocity of the rotating laser beam
	 * @param heat_sigma the laser beam's expansion
	 * @param max_level the max. refinement level
	 * @param refine_threshold threshold for refineing the grid during solution process
	 * @param coarsen_threshold threshold for coarsening the grid during solution process
	 * @param heat the heat factor to initialize the grid
	 */
	LaserHeatEquationSolver2D(double beam_velocity, double heat_sigma, size_t max_level, double refine_threshold, double coarsen_threshold, double heat);

	/**
	 * Std-Destructor of the solver
	 */
	virtual ~LaserHeatEquationSolver2D();

	void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20);

	void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, size_t NumImEul = 0);

	/**
	 * grid initialization for rotating laser test case
	 *
	 * @param alpha reference to the coefficient's vector
	 * @param nRefinements number of initial refinement before solving the heat equation
	 */
	void refineInitialGridWithLaserHeat(sg::base::DataVector& alpha, size_t nRefinements);

	void initScreen();
};

}
}

#endif /* LASERHEATEQUATIONSOLVER2D_HPP */
