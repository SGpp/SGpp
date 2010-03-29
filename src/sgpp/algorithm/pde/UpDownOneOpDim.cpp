/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#include "algorithm/pde/UpDownOneOpDim.hpp"

namespace sg
{

UpDownOneOpDim::UpDownOneOpDim(GridStorage* storage, DataVector& coef)
{
	this->storage = storage;
	this->coefs = &coef;
}

UpDownOneOpDim::UpDownOneOpDim(GridStorage* storage)
{
	this->storage = storage;
	this->coefs = NULL;
}

UpDownOneOpDim::~UpDownOneOpDim()
{
}

void UpDownOneOpDim::mult(DataVector& alpha, DataVector& result)
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
		if (this->coefs != NULL)
		{
			result.axpy_parallel(this->coefs->get(i),beta);
		}
		else
		{
			result.add(beta);
		}
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
			{
				if (this->coefs != NULL)
				{
					result.axpy(this->coefs->get(i),beta);
				}
				else
				{
					result.add(beta);
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
		this->updown(alpha, beta, storage->dim() - 1, i);
		if (this->coefs != NULL)
		{
			result.axpy(this->coefs->get(i),beta);
		}
		else
		{
			result.add(beta);
		}
	}
#endif
}

#ifndef USEOMPTHREE
void UpDownOneOpDim::updown(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim)
{
	if(dim == op_dim)
	{
		specialOP(alpha, result, dim, op_dim);
	}
	else
	{
		//Unidirectional scheme
		if(dim > 0)
		{
			// Reordering ups and downs
			DataVector temp(alpha.getSize());
			up(alpha, temp, dim);
			updown(temp, result, dim-1, op_dim);


			// Same from the other direction:
			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1, op_dim);
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

void UpDownOneOpDim::specialOP(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upOpDim(alpha, temp, dim);
		updown(temp, result, dim-1, op_dim);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, op_dim);
		downOpDim(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		upOpDim(alpha, result, dim);

		DataVector temp(alpha.getSize());
		downOpDim(alpha, temp, dim);

		result.add(temp);
	}
}
#endif

#ifdef USEOMPTHREE
void UpDownOneOpDim::updown_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim)
{
	if(dim == op_dim)
	{
		specialOP_parallel(alpha, result, dim, op_dim);
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
				updown_parallel(temp, result, dim-1, op_dim);
			}

			// Same from the other direction:
			#pragma omp task shared(alpha, temp_two, result_temp)
			{
				updown_parallel(alpha, temp_two, dim-1, op_dim);
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

void UpDownOneOpDim::specialOP_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim)
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
			upOpDim(alpha, temp, dim);
			updown_parallel(temp, result, dim-1, op_dim);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, op_dim);
			downOpDim(temp_two, result_temp, dim);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		upOpDim(alpha, result, dim);

		#pragma omp task shared(alpha, temp)
		downOpDim(alpha, temp, dim);

		#pragma omp taskwait

		result.add(temp);
	}
}
#endif

}
