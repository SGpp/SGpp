/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007 Jörg Blank (blankj@in.tum.de)                          */
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

#include "basis/linear/noboundary/operation/pde/OperationLaplaceLinear.hpp"

#include "basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationLaplaceLinear::OperationLaplaceLinear(GridStorage* storage) : UnidirGradient(storage)
{
}

OperationLaplaceLinear::~OperationLaplaceLinear()
{
}

void OperationLaplaceLinear::mult(DataVector& alpha, DataVector& result)
{
	this->updown(alpha, result);
}

#ifndef USEOMPTHREE
void OperationLaplaceLinear::gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
{
	// In direction gradient_dim we only calculate the norm of the gradient
	// The up-part is empty, thus omitted
	if(dim > 0)
	{
		DataVector temp(alpha.getSize());
		updown(alpha, temp, dim-1, gradient_dim);
		downGradient(temp, result, gradient_dim);
	}
	else
	{
		// Terminates dimension recursion
		downGradient(alpha, result, gradient_dim);
	}
}
#endif

#ifdef USEOMPTHREE
void OperationLaplaceLinear::gradient_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
{
	// In direction gradient_dim we only calculate the norm of the gradient
	// The up-part is empty, thus omitted
	if(dim > 0)
	{
		DataVector temp(alpha.getSize());
		updown_parallel(alpha, temp, dim-1, gradient_dim);
		downGradient(temp, result, gradient_dim);
	}
	else
	{
		// Terminates dimension recursion
		downGradient(alpha, result, gradient_dim);
	}
}
#endif

void OperationLaplaceLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiUpBBLinear func(this->storage);
	sweep<detail::PhiPhiUpBBLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiDownBBLinear func(this->storage);
	sweep<detail::PhiPhiDownBBLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceLinear::downGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// traverse all basis function by sequence number
	for(size_t i = 0; i < storage->size(); i++)
	{
		GridStorage::index_type::level_type level;
		GridStorage::index_type::index_type index;
		(*storage)[i]->get(dim, level, index);
		//only affects the diagonal of the stiffness matrix
		result[i] = alpha[i]*pow(2.0, static_cast<int>(level+1));
	}
}

void OperationLaplaceLinear::upGradient(DataVector& alpha, DataVector& result, size_t dim)
{
}

}
