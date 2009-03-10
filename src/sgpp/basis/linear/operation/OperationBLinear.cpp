/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "basis/basis.hpp"

#include "basis/linear/operation/OperationBLinear.hpp"

#include "sgpp.hpp"

#include "data/DataVector.h"

namespace sg
{
/**
 * Multiplication with vector, not transposed, linear sparse grid
 *
 * @param alpha coefficients of the sparse grid's base functions
 * @param data the vector that should be multiplied
 * @param result the result vector of the matrix vector multiplication
 */
void OperationBLinear::mult(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SLinearBase> op;
	linear_base<unsigned int, unsigned int> base;

	op.mult(storage, base, alpha, data, result);
}

/**
 * Multiplication with vector, transposed, linear sparse grid
 *
 * @param alpha coefficients of the sparse grid's base functions
 * @param data the vector that should be multiplied
 * @param result the result vector of the matrix vector multiplication
 */
void OperationBLinear::multTranspose(DataVector& alpha, DataVector& data, DataVector& result)
{
	AlgorithmB<SLinearBase> op;
	linear_base<unsigned int, unsigned int> base;

	op.mult_transpose(storage, base, alpha, data, result);
}

}
