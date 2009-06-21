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

#ifndef APPLYMATRIX_HPP
#define APPLYMATRIX_HPP

#include "data/DataVector.hpp"
#include "grid/Grid.hpp"

namespace sg
{

/**
 * Virtual class that defines the ApplyMatrix functor, that is used
 * in the Conjugate Gradient solver
 *
 * The derivatives of ApplyMatrix should implement the a matrix functor
 * that does the same as an matrix vector multiplication
 *
 * E.g. if you following equation \f$A\vec{x}=\vec{b}\f$ this functor should
 * modify the vector \f$\vec{x}\f$ in a way like Matrix \f$A\f$ would do.
 */
class ApplyMatrix
{
public:
	/**
	 * Std-Constructor
	 */
	ApplyMatrix() {}

	/**
	 * Std-Destructor
	 */
	virtual ~ApplyMatrix() {}

	/**
	 * Operator, that implements the ApplyMatrix Method in the Conjugate Gradients solver
	 *
	 * @param SparseGrid the grid whose matrix is applied to the vector x
	 * @param x the vector that is multiplied by the matrix
	 * @param b the result of the matrix vector multplication
	 */
	virtual void operator()(DataVector& data, DataVector& x, DataVector& b) = 0;
};

}


#endif /* APPLYMATRIX_HPP */
