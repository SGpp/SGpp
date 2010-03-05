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

#include "basis/linear/noboundary/operation/pde/OperationLTwoDotProductLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/LaplaceDownLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/LaplaceUpLinear.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationLTwoDotProductLinear::OperationLTwoDotProductLinear(GridStorage* storage)
{
	this->storage = storage;
}

OperationLTwoDotProductLinear::~OperationLTwoDotProductLinear()
{
}

void OperationLTwoDotProductLinear::mult(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

	this->updown(alpha, result, storage->dim()-1);
}

/**
 * Recursive procedure for updown
 *
 * @param dim the current dimension
 * @param alpha vector of coefficients
 * @param result vector to store the results in
 */
void OperationLTwoDotProductLinear::updown(DataVector& alpha, DataVector& result, size_t dim)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		up(alpha, temp, dim);
		updown(temp, result, dim-1);

		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1);
		down(temp, result_temp, dim);

		//Overall memory use: 2*|alpha|*(d-1)
		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		up(alpha, result, dim);

		DataVector temp(alpha.getSize());
		down(alpha, temp, dim);

		result.add(temp);
	}
}

void OperationLTwoDotProductLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::LaplaceUpLinear func(this->storage);
	sweep<detail::LaplaceUpLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationLTwoDotProductLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::LaplaceDownLinear func(this->storage);
	sweep<detail::LaplaceDownLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
