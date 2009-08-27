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

#include "basis/lineartrapezoidboundary/operation/finance/OperationDeltaLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiDownLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiUpLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiDownLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiUpLinearTrapezoidBoundary.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationDeltaLinearTrapezoidBoundary::OperationDeltaLinearTrapezoidBoundary(GridStorage* storage, DataVector& mu)
{
	this->storage = storage;
	this->mus = &mu;
}

OperationDeltaLinearTrapezoidBoundary::~OperationDeltaLinearTrapezoidBoundary()
{
}

void OperationDeltaLinearTrapezoidBoundary::mult(DataVector& alpha, DataVector& result)
{
	DataVector beta(result.getSize());
	result.setAll(0.0);

	for(size_t i = 0; i < storage->dim(); i++)
	{
		this->updown(alpha, beta, storage->dim() - 1, i);
		result.axpy((-1.0)*mus->get(i),beta);
	}
}

void OperationDeltaLinearTrapezoidBoundary::updown(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
{
	if(dim == gradient_dim)
	{
		gradient(alpha, result, dim, gradient_dim);
	}
	else
	{
		//Unidirectional scheme
		if(dim > 0)
		{
			// Reordering ups and downs
			// Use previously calculated ups for all future calculations
			// U* -> UU* and UD*

			DataVector temp(alpha.getSize());
			up(alpha, temp, dim);
			updown(temp, result, dim-1, gradient_dim);


			// Same from the other direction:
			// *D -> *UD and *DD

			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1, gradient_dim);
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
}

void OperationDeltaLinearTrapezoidBoundary::gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upGradient(alpha, temp, dim);
		updown(temp, result, dim-1, gradient_dim);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, gradient_dim);
		downGradient(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		upGradient(alpha, result, dim);

		DataVector temp(alpha.getSize());
		downGradient(alpha, temp, dim);

		result.add(temp);
	}
}

void OperationDeltaLinearTrapezoidBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiUpLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearTrapezoidBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiDownLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearTrapezoidBoundary::upGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiUpLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XdPhiPhiUpLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearTrapezoidBoundary::downGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiDownLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XdPhiPhiDownLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
