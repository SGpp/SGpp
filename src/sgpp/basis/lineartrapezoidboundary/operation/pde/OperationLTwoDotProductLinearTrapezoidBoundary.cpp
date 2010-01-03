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

#include "basis/lineartrapezoidboundary/operation/pde/OperationLTwoDotProductLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "algorithm/common/sweep.hpp"

namespace sg
{

OperationLTwoDotProductLinearTrapezoidBoundary::OperationLTwoDotProductLinearTrapezoidBoundary(GridStorage* storage)
{
	this->storage = storage;
}

OperationLTwoDotProductLinearTrapezoidBoundary::~OperationLTwoDotProductLinearTrapezoidBoundary()
{
}

void OperationLTwoDotProductLinearTrapezoidBoundary::mult(DataVector& alpha, DataVector& result)
{
	DataVector beta(result.getSize());
	result.setAll(0.0);
#ifdef USEOMP
#ifdef USEOMPTHREE
	#pragma omp parallel
	{
#ifndef AIX_XLC
		#pragma omp single nowait
#endif
#ifdef AIX_XLC
		#pragma omp single
#endif
		{
			this->updown_parallel(alpha, beta, storage->dim() - 1);
		}
	}
	result.add(beta);
#endif
#ifndef USEOMPTHREE
	this->updown(alpha, beta, storage->dim() - 1);

	result.add(beta);
#endif
#endif
#ifndef USEOMP
	this->updown(alpha, beta, storage->dim() - 1);

	result.add(beta);
#endif
}

#ifndef USEOMPTHREE
void OperationLTwoDotProductLinearTrapezoidBoundary::updown(DataVector& alpha, DataVector& result, size_t dim)
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
#endif

#ifdef USEOMPTHREE
void OperationLTwoDotProductLinearTrapezoidBoundary::updown_parallel(DataVector& alpha, DataVector& result, size_t dim)
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
			updown_parallel(temp, result, dim-1);
		}


		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1);
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
#endif

void OperationLTwoDotProductLinearTrapezoidBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLTwoDotProductLinearTrapezoidBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	// phi * phi
	detail::PhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

	s.sweep1D_Boundary(alpha, result, dim);
}

}
