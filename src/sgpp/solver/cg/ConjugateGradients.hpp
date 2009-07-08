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

#ifndef CONJUGATEGRADIENTS_HPP
#define CONJUGATEGRADIENTS_HPP

#include "algorithm/classification/DMSystemMatrix.hpp"
#include "solver/SGSolver.hpp"
#include "data/DataVector.hpp"

namespace sg
{

class ConjugateGradients : public SGSolver
{
private:


public:
	/**
	 * Std-Constructor
	 */
	ConjugateGradients(size_t imax, double epsilon);

	/**
	 * Std-Destructor
	 */
	virtual ~ConjugateGradients();

	virtual void solve(OperationMatrix& SystemMatrix, DataVector& alpha, DataVector& b, bool reuse = false, bool verbose = false, double max_threshold = -1.0);

	// Define functions for observer pattern in python

	/**
	 * function that signals the start of the CG method (used in python)
	 */
	virtual void starting();

	/**
	 * function that signals the start of the calculation of the CG method (used in python)
	 */
	virtual void calcStarting();

	/**
	 * function that signals that one iteration step of the CG method has been completed (used in python)
	 */
	virtual void iterationComplete();

	/**
	 * function that signals the finish of the cg method (used in python)
	 */
	virtual void complete();
};

}

#endif /* CONJUGATEGRADIENTS_HPP */
