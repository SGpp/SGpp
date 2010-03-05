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

#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationDeltaLinearBoundary::OperationDeltaLinearBoundary(GridStorage* storage, DataVector& coef)
{
	this->storage = storage;
	this->coefs = &coef;
}

OperationDeltaLinearBoundary::~OperationDeltaLinearBoundary()
{
}

void OperationDeltaLinearBoundary::mult(DataVector& alpha, DataVector& result)
{
	result.setAll(0.0);

#ifdef USEOMP
#ifdef USEOMPTHREE
	DataVector beta(result.getSize());

	for(size_t i = 0; i < storage->dim(); i++)
	{
		#pragma omp parallel
		{
#ifndef AIX_XLC
			#pragma omp single nowait
#endif
#ifdef AIX_XLC
			#pragma omp single
#endif
			{
				this->updown_parallel(alpha, beta, storage->dim() - 1, i);
			}
		}
		result.axpy_parallel(this->coefs->get(i),beta);
	}
#endif
#ifndef USEOMPTHREE
	#pragma omp parallel shared(result)
	{
		#pragma omp for schedule(static)
		for(size_t i = 0; i < storage->dim(); i++)
		{
			DataVector beta(result.getSize());

			this->updown(alpha, beta, storage->dim() - 1, i);

			#pragma omp critical
			result.axpy(this->coefs->get(i),beta);
		}
	}
#endif
#endif
#ifndef USEOMP
	DataVector beta(result.getSize());

	for(size_t i = 0; i < storage->dim(); i++)
	{
		this->updown(alpha, beta, storage->dim() - 1, i);
		result.axpy(this->coefs->get(i),beta);
	}
#endif
}

#ifndef USEOMPTHREE
void OperationDeltaLinearBoundary::updown(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
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
			DataVector temp(alpha.getSize());
			up(alpha, temp, dim);
			updown(temp, result, dim-1, gradient_dim);


			// Same from the other direction:
			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1, gradient_dim);
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

void OperationDeltaLinearBoundary::gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
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
#endif

#ifdef USEOMPTHREE
void OperationDeltaLinearBoundary::updown_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
{
	if(dim == gradient_dim)
	{
		gradient_parallel(alpha, result, dim, gradient_dim);
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
				updown_parallel(temp, result, dim-1, gradient_dim);
			}

			// Same from the other direction:
			#pragma omp task shared(alpha, temp_two, result_temp)
			{
				updown_parallel(alpha, temp_two, dim-1, gradient_dim);
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

void OperationDeltaLinearBoundary::gradient_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
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
			updown_parallel(temp, result, dim-1, gradient_dim);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, gradient_dim);
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
#endif

void OperationDeltaLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::upGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::XPhidPhiUpBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationDeltaLinearBoundary::downGradient(DataVector& alpha, DataVector& result, size_t dim)
{
	// x * phi * dphi
	detail::XPhidPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::XPhidPhiDownBBLinearBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
