/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#ifndef EXPLICITEULER_HPP
#define EXPLICITEULER_HPP

#include "operation/common/OperationMatrix.hpp"
#include "solver/ODESolver.hpp"

namespace sg
{

/**
 * This class implements the explicit Euler method
 * for solving ordinary partial equations
 */
class ExplicitEuler : public ODESolver
{
private:
	/// the number of CG maximum CG iterations
	size_t maxCGIterations;
	/// the CG's epsilon
	double epsilonCG;
	/// specifies if a grid evaluation should be execute in every time step
	bool bAnimation;

public:
	/**
	 * Std-Constructer
	 *
	 * @param imax number of maximum executed iterations
	 * @param timestepSize the size of one timestep
	 * @param generateAnimation set this, if you want to create a grid evaluation in every time step, in order to create an animation
	 */
	ExplicitEuler(size_t imax, double timestepSize, size_t iMaxCG, double epsilonCG, bool generateAnimation = false);

	/**
	 * Std-Destructor
	 */
	virtual ~ExplicitEuler();

	virtual void solve(OperationSolverMatrix& SystemMatrix, DataVector& alpha, bool verbose = false);
};

}

#endif /* EXPLICITEULER_HPP */
