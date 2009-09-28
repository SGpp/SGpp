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

#include "basis/lineartrapezoidboundary/operation/pde/finance/OperationGammaPartThreeLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/SqXdPhidPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/SqXdPhidPhiUpBBLinearTrapezoidBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>

namespace sg
{

OperationGammaPartThreeLinearTrapezoidBoundary::OperationGammaPartThreeLinearTrapezoidBoundary(GridStorage* storage, DataVector& sigma, DataVector& rho)
{
	this->storage = storage;
	this->sigmas = &sigma;
	this->rhos = &rho;
}

OperationGammaPartThreeLinearTrapezoidBoundary::~OperationGammaPartThreeLinearTrapezoidBoundary()
{
}

void OperationGammaPartThreeLinearTrapezoidBoundary::mult(DataVector& alpha, DataVector& result)
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
				//std::cout << beta.toString() << std::endl;
				result.axpy((0.5)*sigmas->get(i)*sigmas->get(j)*rhos->get((storage->dim()*i)+j),beta);
			}
			// Use the symmetry of the operation (i,j)+(j,i) = 2*(i,j)
			if (j < i)
			{
				this->updown(alpha, beta, storage->dim() - 1, i, j);
				result.axpy((1.0)*sigmas->get(i)*sigmas->get(j)*rhos->get((storage->dim()*i)+j),beta);
			}
			//std::cout << sigmas->get(i) << " " << sigmas->get(j) << " " << rhos->get((storage->dim()*i)+j) << std::endl;
			//std::cout << result.toString() << std::endl;
		}
	}
}

void OperationGammaPartThreeLinearTrapezoidBoundary::updown(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	if((dim == gradient_dim_one) && (dim == gradient_dim_two))
	{
		gradientSquared(alpha, result, dim, gradient_dim_one, gradient_dim_two);
	}
	else if ((dim == gradient_dim_one || dim == gradient_dim_two) && (gradient_dim_one != gradient_dim_two))
	{
		gradient(alpha, result, dim, gradient_dim_one, gradient_dim_two);
	}
	else
	{
		//Unidirectional scheme
		if(dim > 0)
		{
			// Reordering ups and downs
			DataVector temp(alpha.getSize());
			up(alpha, temp, dim);
			updown(temp, result, dim-1, gradient_dim_one, gradient_dim_two);


			// Same from the other direction:
			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1, gradient_dim_one, gradient_dim_two);
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

void OperationGammaPartThreeLinearTrapezoidBoundary::gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upGradient(alpha, temp, dim);
		updown(temp, result, dim-1, gradient_dim_one, gradient_dim_two);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, gradient_dim_one, gradient_dim_two);
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

void OperationGammaPartThreeLinearTrapezoidBoundary::gradientSquared(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upGradientSquared(alpha, temp, dim);
		updown(temp, result, dim-1, gradient_dim_one, gradient_dim_two);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, gradient_dim_one, gradient_dim_two);
		downGradientSquared(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		upGradientSquared(alpha, result, dim);

		DataVector temp(alpha.getSize());
		downGradientSquared(alpha, temp, dim);

		result.add(temp);
	}
}

void OperationGammaPartThreeLinearTrapezoidBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartThreeLinearTrapezoidBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartThreeLinearTrapezoidBoundary::upGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XdPhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartThreeLinearTrapezoidBoundary::downGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XdPhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartThreeLinearTrapezoidBoundary::upGradientSquared(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	detail::SqXdPhidPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::SqXdPhidPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaPartThreeLinearTrapezoidBoundary::downGradientSquared(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	detail::SqXdPhidPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::SqXdPhidPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
