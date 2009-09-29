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

#ifndef DIRICHLETUPDATEVECTOR_HPP
#define DIRICHLETUPDATEVECTOR_HPP

#include "grid/common/BoundingBox.hpp"
#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

/**
 * This class is useful if you do some PDE calculations with Dirichlet Boundary
 * Conditions. Doing this, e.g. you might wish to add some solution from a timestep to
 * the current coefficients of the sparse grid. If you are using Dirichlet
 * conditions you mustn't overwrite the values on the boundaries in your coefficient
 * vector.
 *
 * This class implements a functor that uses the Bounding Box of the grid to determine, if
 * a boundary has to implement Dirichlet boundary conditions. In that case, simply use this
 * to replace all values in the update vector on these boundaries with zero, so you can
 * safely add the resulting vector to your solution.
 */
class DirichletUpdateVector
{
private:
	/// bounding box of the grid
	BoundingBox* myBoundingBox;
	/// Grid Storage object
	GridStorage* storage;

public:
	/**
	 * Std-Constructor
	 *
	 * @param storage the grid's storage object; needed to determine the bounding box and to iterate of the entries in the coefficient vector
	 */
	DirichletUpdateVector(GridStorage* storage);

	/**
	 * Std-Destructor
	 */
	~DirichletUpdateVector();

	/**
	 * Replace the boundary entries in updateVector with the one from sourceVector
	 * only in that dimension, for which Dirichlet Boundary Conditions
	 * were specified
	 *
	 * @param updateVector the vector that should be updated
	 * @param sourceVector the vector that contains the correct boundary values
	 */
	void applyDirichletConditions(DataVector& updateVector, DataVector& sourceVector);

	/**
	 * Replace the boundary entries in updateVector with Zero only in that dimension, for which Dirichlet Boundary Conditions
	 * were specified
	 *
	 * @param updateVector the vector that should be updated
	 */
	void setBoundariesToZero(DataVector& updateVector);

	/**
	 * Multiplies the values on the boundary with a constant value
	 *
	 * @param updateVector the vector that should be updated
	 * @param value the value that is multiplied with the value on the boundaries
	 */
	void multiplyBoundary(DataVector& updateVector, double value);
};

}

#endif /* DIRICHLETUPDATEVECTOR_HPP */
