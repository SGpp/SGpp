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

#ifndef SGSOLVER_HPP
#define SGSOLVER_HPP

namespace sg
{

/**
 * Abstract class that defines a solver used in Sparse Grids
 * Applications
 */
class SGSolver
{
private:
	/// Number of Iterations needed for the solve
	size_t nIterations;
	/// Number of maximum iterations for cg
	size_t nMaxIterations;
	/// final residuum
	double finalResiduum;

public:
	/**
	 * Std-Constructor
	 */
	SGSolver(size_t nMaximumIterations) : nMaxIterations(nMaximumIterations), nIterations(0), finalResiduum(0.0) { }

	/**
	 * Std-Destructor
	 */
	virtual ~SGSolver() { }

	/**
	 * Pure virtual Function that defines a solve method
	 */
	virtual void solve() = 0;

	/**
	 * function that returns the number of needed solve steps
	 *
	 * @return the number of needed solve steps of the sovler
	 */
	size_t getNumberIterations()
	{
		return nIterations;
	}

	/**
	 * function the returns the final residuum, error of the solver
	 *
	 * @sreturn the final residuum
	 */
	double getFinalResiduum()
	{
		return finalResiduum;
	}
};

}

#endif /* SGSOLVER_HPP */
