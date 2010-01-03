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

#ifndef EULER_HPP
#define EULER_HPP

#include "application/common/ScreenOutput.hpp"
#include "solver/ODESolver.hpp"
#include <string>

namespace sg
{

/**
 * This class implements the explicit Euler method
 * for solving ordinary partial equations
 *
 * @version $HEAD$
 */
class Euler : public ODESolver
{
private:
	/// the number of CG maximum CG iterations
	size_t maxCGIterations;
	/// the CG's epsilon
	double epsilonCG;
	/// specifies if a grid evaluation should be execute in every time step
	bool bAnimation;
	/// specifies the evaluation per dimension when a animation is created
	size_t evalsAnimation;
	/// specifies the type of euler that should be executed
	std::string ExMode;
	/// Pointer to ScreenOutput object
	ScreenOutput* myScreen;

public:
	/**
	 * Std-Constructer
	 *
	 * @param Mode the mode of the euler that should be executed, must be ExEul or ImEul
	 * @param imax number of maximum executed iterations
	 * @param timestepSize the size of one timestep
	 * @param iMaxCG maximum number of CG steps
	 * @param epsilonCG the epsilon used in CG
	 * @param generateAnimation set this, if you want to create a grid evaluation in every time step, in order to create an animation
	 * @param numEvalsAnimation specifies the evaluation per dimension when a animation is created
	 * @param screen possible pointer to a ScreenOutput object
	 */
	Euler(std::string Mode, size_t imax, double timestepSize, size_t iMaxCG, double epsilonCG, bool generateAnimation = false, size_t numEvalsAnimation = 20, ScreenOutput* screen = NULL);

	/**
	 * Std-Destructor
	 */
	virtual ~Euler();

	virtual void solve(OperationODESolverMatrix& SystemMatrix, DataVector& alpha, bool verbose = false);
};

}

#endif /* EULER_HPP */
