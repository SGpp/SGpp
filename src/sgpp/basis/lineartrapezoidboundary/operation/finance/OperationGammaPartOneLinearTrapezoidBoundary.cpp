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

#include "basis/lineartrapezoidboundary/operation/finance/OperationGammaPartOneLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/XSurfaceIntegralBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/SqXSurfaceIntegralBBLinearTrapezoidBoundary.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationGammaPartOneLinearTrapezoidBoundary::OperationGammaPartOneLinearTrapezoidBoundary(GridStorage* storage, DataVector& sigma, DataVector& rho)
{
	this->storage = storage;
	this->sigmas = &sigma;
	this->rhos = &rho;
}

OperationGammaPartOneLinearTrapezoidBoundary::~OperationGammaPartOneLinearTrapezoidBoundary()
{
}

void OperationGammaPartOneLinearTrapezoidBoundary::mult(DataVector& alpha, DataVector& result)
{
	DataVector beta(result.getSize());
	result.setAll(0.0);

	for(size_t i = 0; i < storage->dim(); i++)
	{
		for(size_t j = 0; j < storage->dim(); j++)
		{
			// Calculate the "diagonal" of the operation
			if (j == i)
			{
				this->updown(alpha, beta, storage->dim() - 1, i, j);
				result.axpy((0.5)*sigmas->get(i)*sigmas->get(j)*rhos->get((storage->dim()*i)+j),beta);
			}
			// Use the symmetry of the operation (i,j)+(j,i) = 2*(i,j)
			if (j < i)
			{
				this->updown(alpha, beta, storage->dim() - 1, i, j);
				result.axpy((1.0)*sigmas->get(i)*sigmas->get(j)*rhos->get((storage->dim()*i)+j),beta);
			}
		}
	}
}

void OperationGammaPartOneLinearTrapezoidBoundary::updown(DataVector& alpha, DataVector& result, size_t dim, size_t operation_dim_one, size_t operation_dim_two)
{
	if((dim == operation_dim_one) && (dim == operation_dim_two))
	{
		SurfaceIntegralSquared(alpha, result, dim, operation_dim_one, operation_dim_two);
	}
	else if (dim == operation_dim_two)
	{
		SurfaceIntegral(alpha, result, dim, operation_dim_one, operation_dim_two);
	}
	else if (dim == operation_dim_one)
	{
		gradient(alpha, result, dim, operation_dim_one, operation_dim_two);
	}
	else
	{
		//Unidirectional scheme
		if(dim > 0)
		{
			// Reordering ups and downs
			DataVector temp(alpha.getSize());
			up(alpha, temp, dim);
			updown(temp, result, dim-1, operation_dim_one, operation_dim_two);

			// Same from the other direction:
			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1, operation_dim_one, operation_dim_two);
			down(temp, result_temp, dim);

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

void OperationGammaPartOneLinearTrapezoidBoundary::gradient(DataVector& alpha, DataVector& result, size_t dim, size_t operation_dim_one, size_t operation_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upGradient(alpha, temp, dim);
		updown(temp, result, dim-1, operation_dim_one, operation_dim_two);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, operation_dim_one, operation_dim_two);
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

void OperationGammaPartOneLinearTrapezoidBoundary::SurfaceIntegral(DataVector& alpha, DataVector& result, size_t dim, size_t operation_dim_one, size_t operation_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		updown(temp, result, dim-1, operation_dim_one, operation_dim_two);

		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, operation_dim_one, operation_dim_two);
		calcSurfaceIntegral(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		calcSurfaceIntegral(alpha, result, dim);
	}
}

void OperationGammaPartOneLinearTrapezoidBoundary::SurfaceIntegralSquared(DataVector& alpha, DataVector& result, size_t dim, size_t operation_dim_one, size_t operation_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		updown(temp, result, dim-1, operation_dim_one, operation_dim_two);

		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, operation_dim_one, operation_dim_two);
		calcSurfaceIntegralSquared(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		calcSurfaceIntegralSquared(alpha, result, dim);
	}

}

void OperationGammaPartOneLinearTrapezoidBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartOneLinearTrapezoidBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartOneLinearTrapezoidBoundary::upGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XdPhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartOneLinearTrapezoidBoundary::downGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XdPhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartOneLinearTrapezoidBoundary::calcSurfaceIntegral(DataVector& alpha, DataVector& result, size_t dim)
{
	// Surface Integral
	detail::XSurfaceIntegralBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XSurfaceIntegralBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartOneLinearTrapezoidBoundary::calcSurfaceIntegralSquared(DataVector& alpha, DataVector& result, size_t dim)
{
	// Surface Integral
	detail::SqXSurfaceIntegralBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::SqXSurfaceIntegralBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
