/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef POISSONEQUATIONSOLVERMPI_HPP
#define POISSONEQUATIONSOLVERMPI_HPP

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
using namespace sg::base;

namespace sg
{

/**
 * This class provides a simple-to-use solver of the multi dimensional
 * Poisson Equation on Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the
 * Poisson Equation on Sparse Grids!
 *
 * This version offers support for MPI!
 *
 * @version $HEAD$
 */
class PoissonEquationSolverMPI : public EllipticPDESolver
{
private:
	/// screen object used in this solver
	ScreenOutput* myScreen;

public:
	/**
	 * Std-Constructor of the solver
	 */
	PoissonEquationSolverMPI();

	/**
	 * Std-Destructor of the solver
	 */
	virtual ~PoissonEquationSolverMPI();

	void constructGrid(BoundingBox& myBoundingBox, size_t level);

	void solvePDE(DataVector& alpha, DataVector& rhs, size_t maxCGIterations, double epsilonCG, bool verbose = false);

	/**
	 * Inits the grid with a smooth heat distribution (based on
	 * a std-normal distribution) on its boundaries
	 *
	 * @param alpha reference to the coefficients vector
	 * @param mu the exspected value of the normal distribution
	 * @param sigma the sigma of the normal distribution
	 * @param factor a factor that is used to stretch the function values
	 */
	void initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma, double factor);

	/**
	 * Inits the screen object
	 */
	void initScreen();
};

}

#endif /* POISSONEQUATIONSOLVERMPI_HPP */
