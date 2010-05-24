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

#ifndef PARABOLICPDESOLVER_HPP
#define PARABOLICPDESOLVER_HPP

#include "application/pde/PDESolver.hpp"

namespace sg
{

/**
 * This class extends the PDESolver with functions that are needed to
 * solve parabolic PDEs
 *
 * @version $HEAD$
 */
class ParabolicPDESolver : public PDESolver
{
protected:
	/// the size of one timestep
	//double timestepSize;
	/// The number of timesteps that are executed during solving
	//size_t nTimesteps;

public:
	/**
	 * Std-Constructor of the solver
	 */
	ParabolicPDESolver();

	/**
	 * Std-Destructor of the solver
	 */
	virtual ~ParabolicPDESolver();

	/**
	 * Call this routine to use an explicit Euler algorithm to solve the parabolic PDE
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
	virtual void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20) = 0;

	/**
	 * Call this routine to use an explicit Euler algorithm to solve the parabolic PDE
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
	virtual void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose = false, bool generateAnimation = false, size_t numEvalsAnimation = 20) = 0;

	/**
	 * Call this routine to use the Crank Nicolson algorithm to solve the parabolic PDE
	 *
	 * @param numTimesteps the number of timesteps that should be executed
	 * @param timestepsize the size of the interval one timestep moves forward
	 * @param maxCGIterations the maximum of interation in the CG solver
	 * @param epsilonCG the epsilon used in the CG
	 * @param alpha the coefficients of the Sparse Gird's basis functions
	 * @param NumImEul specifies how many ImEul steps should be executed before CrNic is used, default is 3
	 */
	virtual void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul = 3) = 0;
};

}

#endif /* PARABOLICPDESOLVER_HPP */
