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

#include "basis/lineartrapezoidboundary/operation/pde/finance/OperationGammaLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/XPhidPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/XPhidPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/SqXdPhidPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/SqXdPhidPhiUpBBLinearTrapezoidBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include <iostream>

namespace sg
{

OperationGammaLinearTrapezoidBoundary::OperationGammaLinearTrapezoidBoundary(GridStorage* storage, DataVector& coef)
{
	this->storage = storage;
	this->coefs = &coef;
}

OperationGammaLinearTrapezoidBoundary::~OperationGammaLinearTrapezoidBoundary()
{
}

void OperationGammaLinearTrapezoidBoundary::mult(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);
#ifdef USEOMP
#ifdef USEOMPTHREE
	DataVector beta(result.getSize());

	for(size_t i = 0; i < storage->dim(); i++)
	{
		for(size_t j = 0; j < storage->dim(); j++)
		{
			// use the operator's symmetry
			if ( j <= i)
			{
				#pragma omp parallel
				{
					#pragma omp single nowait
					{
						this->updown_parallel(alpha, beta, storage->dim() - 1, i, j);
					}
				}
			}
			// Calculate the "diagonal" of the operation
			if (j == i)
			{
				result.axpy_parallel(this->coefs->get((storage->dim()*i)+j),beta);
			}
			// Use the symmetry of the operation (i,j)+(j,i) = 2*(i,j)
			if (j < i)
			{
				result.axpy_parallel(this->coefs->get((storage->dim()*i)+j),beta);
			}
		}
	}
#endif
#ifndef USEOMPTHREE
	#pragma omp parallel shared(result)
	{
		#pragma omp for schedule(static)
		for(size_t i = 0; i < storage->dim(); i++)
		{
			for(size_t j = 0; j < storage->dim(); j++)
			{
				DataVector beta(result.getSize());

				// use the operator's symmetry
				if ( j <= i)
				{
					this->updown(alpha, beta, storage->dim() - 1, i, j);
				}

				// Calculate the "diagonal" of the operation
				if (j == i)
				{
					#pragma omp critical
					result.axpy(this->coefs->get((storage->dim()*i)+j),beta);
				}
				// Use the symmetry of the operation (i,j)+(j,i) = 2*(i,j)
				if (j < i)
				{
					#pragma omp critical
					result.axpy(this->coefs->get((storage->dim()*i)+j),beta);
				}
			}
		}
	}
#endif
#endif
#ifndef USEOMP
	DataVector beta(result.getSize());

	for(size_t i = 0; i < storage->dim(); i++)
	{
		for(size_t j = 0; j < storage->dim(); j++)
		{
			// use the operator's symmetry
			if ( j <= i)
			{
				this->updown(alpha, beta, storage->dim() - 1, i, j);
			}
			// Calculate the "diagonal" of the operation
			if (j == i)
			{
				result.axpy(this->coefs->get((storage->dim()*i)+j),beta);
			}
			// Use the symmetry of the operation (i,j)+(j,i) = 2*(i,j)
			if (j < i)
			{
				result.axpy(this->coefs->get((storage->dim()*i)+j),beta);
			}
		}
	}
#endif
}

#ifndef USEOMPTHREE
void OperationGammaLinearTrapezoidBoundary::updown(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	if((dim == gradient_dim_one) && (dim == gradient_dim_two))
	{
		gradientSquared(alpha, result, dim, gradient_dim_one, gradient_dim_two);
	}
	else if ((dim == gradient_dim_one || dim == gradient_dim_two) && (gradient_dim_one != gradient_dim_two))
	{
		if (dim == gradient_dim_one)
		{
			gradient(alpha, result, dim, gradient_dim_one, gradient_dim_two);
		}
		if (dim == gradient_dim_two)
		{
			gradientTestFct(alpha, result, dim, gradient_dim_one, gradient_dim_two);
		}
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

void OperationGammaLinearTrapezoidBoundary::gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
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

void OperationGammaLinearTrapezoidBoundary::gradientTestFct(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upGradientTestFct(alpha, temp, dim);
		updown(temp, result, dim-1, gradient_dim_one, gradient_dim_two);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, gradient_dim_one, gradient_dim_two);
		downGradientTestFct(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		upGradientTestFct(alpha, result, dim);

		DataVector temp(alpha.getSize());
		downGradientTestFct(alpha, temp, dim);

		result.add(temp);
	}
}

void OperationGammaLinearTrapezoidBoundary::gradientSquared(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
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
#endif

#ifdef USEOMPTHREE
void OperationGammaLinearTrapezoidBoundary::updown_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	if((dim == gradient_dim_one) && (dim == gradient_dim_two))
	{
		gradientSquared_parallel(alpha, result, dim, gradient_dim_one, gradient_dim_two);
	}
	else if ((dim == gradient_dim_one || dim == gradient_dim_two) && (gradient_dim_one != gradient_dim_two))
	{
		gradient_parallel(alpha, result, dim, gradient_dim_one, gradient_dim_two);
	}
	else
	{
		//Unidirectional scheme
		if(dim > 0)
		{
			// Reordering ups and downs
			DataVector temp(alpha.getSize());
			DataVector result_temp(alpha.getSize());
			DataVector temp_two(alpha.getSize());

			#pragma omp task shared(alpha, temp, result)
			{
				up(alpha, temp, dim);
				updown_parallel(temp, result, dim-1, gradient_dim_one, gradient_dim_two);
			}

			// Same from the other direction:
			#pragma omp task shared(alpha, temp_two, result_temp)
			{
				updown_parallel(alpha, temp_two, dim-1, gradient_dim_one, gradient_dim_two);
				down(temp_two, result_temp, dim);
			}

			#pragma omp taskwait

			result.add(result_temp);
		}
		else
		{
			// Terminates dimension recursion
			DataVector temp(alpha.getSize());

			#pragma omp task shared(alpha, result)
			up(alpha, result, dim);

			#pragma omp task shared(alpha, temp)
			down(alpha, temp, dim);

			#pragma omp taskwait

			result.add(temp);
		}
	}
}

void OperationGammaLinearTrapezoidBoundary::gradient_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		DataVector result_temp(alpha.getSize());
		DataVector temp_two(alpha.getSize());

		#pragma omp task shared(alpha, temp, result)
		{
			upGradient(alpha, temp, dim);
			updown_parallel(temp, result, dim-1, gradient_dim_one, gradient_dim_two);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, gradient_dim_one, gradient_dim_two);
			downGradient(temp_two, result_temp, dim);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		upGradient(alpha, result, dim);

		#pragma omp task shared(alpha, temp)
		downGradient(alpha, temp, dim);

		#pragma omp taskwait

		result.add(temp);
	}
}

void OperationGammaLinearTrapezoidBoundary::gradientTestFct_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		DataVector result_temp(alpha.getSize());
		DataVector temp_two(alpha.getSize());

		#pragma omp task shared(alpha, temp, result)
		{
			upGradientTestFct(alpha, temp, dim);
			updown_parallel(temp, result, dim-1, gradient_dim_one, gradient_dim_two);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, gradient_dim_one, gradient_dim_two);
			downGradientTestFct(temp_two, result_temp, dim);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		upGradientTestFct(alpha, result, dim);

		#pragma omp task shared(alpha, temp)
		downGradientTestFct(alpha, temp, dim);

		#pragma omp taskwait

		result.add(temp);
	}
}

void OperationGammaLinearTrapezoidBoundary::gradientSquared_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim_one, size_t gradient_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		DataVector result_temp(alpha.getSize());
		DataVector temp_two(alpha.getSize());

		#pragma omp task shared(alpha, temp, result)
		{
			upGradientSquared(alpha, temp, dim);
			updown_parallel(temp, result, dim-1, gradient_dim_one, gradient_dim_two);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, gradient_dim_one, gradient_dim_two);
			downGradientSquared(temp_two, result_temp, dim);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		upGradientSquared(alpha, result, dim);

		#pragma omp task shared(alpha, temp)
		downGradientSquared(alpha, temp, dim);

		#pragma omp taskwait

		result.add(temp);
	}
}
#endif

void OperationGammaLinearTrapezoidBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearTrapezoidBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearTrapezoidBoundary::upGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XPhidPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearTrapezoidBoundary::downGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XPhidPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearTrapezoidBoundary::upGradientTestFct(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XdPhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearTrapezoidBoundary::downGradientTestFct(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * dphi * phi
	detail::XdPhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::XdPhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearTrapezoidBoundary::upGradientSquared(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	detail::SqXdPhidPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::SqXdPhidPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLinearTrapezoidBoundary::downGradientSquared(DataVector& alpha, DataVector& result, size_t dim)
{
	// x^2 * dphi * dphi
	detail::SqXdPhidPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::SqXdPhidPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
