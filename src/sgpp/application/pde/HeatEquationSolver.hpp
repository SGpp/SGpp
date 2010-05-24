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
 * Heat Equation that uses Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the
 * Heat Equation with Sparse Grids!
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

	void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul = 3);

	/**
	 * This method sets the heat coefficient of the regarded material
	 *
	 * @param a the heat coefficient
	 */
	void setHeatCoefficient(double a);

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
