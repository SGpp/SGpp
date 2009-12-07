/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
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

#ifndef OPERATIONODESOLVERMATRIX_HPP
#define OPERATIONODESOLVERMATRIX_HPP

#include "grid/Grid.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"

namespace sg
{

/**
 * Defines a Systemmatrix that is used to solve parabolic partial
 * differential equations. So an instance of this class has to pass to
 * any ODE Solver used in SGpp.
 *
 * This class defines an elliptic problem in every timestep which is solved
 * using an iterative SLE solver, that solving step is integrated in the
 * ODE Solver.
 */
class OperationODESolverMatrix : public OperationMatrix
{
public:
	/**
	 * Constructor
	 */
	OperationODESolverMatrix() {}

	/**
	 * Destructor
	 */
	virtual ~OperationODESolverMatrix() {}

	/**
	 * Multiplicates a vector with the matrix
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 * @param result DataVector into which the result of the Laplace operation is stored
	 */
	virtual void mult(DataVector& alpha, DataVector& result) = 0;

	/**
	 * generates the right hand side of the system
	 *
	 * @param alpha data that is needed to generate the RHS of the system of linear equations
	 * @param rhs DataVector the contains the RHS
	 */
	virtual void generateRHS(DataVector& alpha, DataVector& rhs) = 0;

	/**
	 * performs some action that might be needed after a timestep has be finished in the ODE
	 * Solver, e.g. some boundary adjustments.
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 */
	virtual void finishTimestep(DataVector& alpha) = 0;

	/**
	 * performs some action that might be needed before a timestep is executed in the ODE
	 * Solver, e.g. some boundary adjustments.
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 */
	virtual void startTimestep(DataVector& alpha) = 0;

	/**
	 * get the pointer to the underlying grid object
	 *
	 * @return returns a pointer to the underlying grid object
	 */
	virtual Grid* getGrid() = 0;
};

}

#endif /* OPERATIONODESOLVERMATRIX_HPP */
