/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/UpDownOneOpDim.hpp"

namespace sg
{

UpDownOneOpDim::UpDownOneOpDim(GridStorage* storage, DataVector& coef)
{
	this->storage = storage;
	this->coefs = &coef;
	this->algoDims = this->storage->getAlgorithmicDimensions();
}

UpDownOneOpDim::UpDownOneOpDim(GridStorage* storage)
{
	this->storage = storage;
	this->coefs = NULL;
	this->algoDims = this->storage->getAlgorithmicDimensions();
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

	for(size_t i = 0; i < this->algoDims.size(); i++)
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
				this->updown_parallel(alpha, beta, this->algoDims.size() - 1, i);
			}
		}
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
#ifndef USEOMPTHREE
	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for(size_t i = 0; i < this->algoDims.size(); i++)
		{
			DataVector beta(result.getSize());

			this->updown(alpha, beta, this->algoDims.size() - 1, i);

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

	for(size_t i = 0; i < this->algoDims.size(); i++)
	{
		this->updown(alpha, beta, this->algoDims.size() - 1, i);
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
			up(alpha, temp, this->algoDims[dim]);
			updown(temp, result, dim-1, op_dim);


			// Same from the other direction:
			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1, op_dim);
			down(temp, result_temp, this->algoDims[dim]);

			result.add(result_temp);
		}
		else
		{
			// Terminates dimension recursion
			up(alpha, result, this->algoDims[dim]);

			DataVector temp(alpha.getSize());
			down(alpha, temp, this->algoDims[dim]);

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
		upOpDim(alpha, temp, this->algoDims[dim]);
		updown(temp, result, dim-1, op_dim);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, op_dim);
		downOpDim(temp, result_temp, this->algoDims[dim]);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		upOpDim(alpha, result, this->algoDims[dim]);

		DataVector temp(alpha.getSize());
		downOpDim(alpha, temp, this->algoDims[dim]);

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
				up(alpha, temp, this->algoDims[dim]);
				updown_parallel(temp, result, dim-1, op_dim);
			}

			// Same from the other direction:
			#pragma omp task shared(alpha, temp_two, result_temp)
			{
				updown_parallel(alpha, temp_two, dim-1, op_dim);
				down(temp_two, result_temp, this->algoDims[dim]);
			}

			#pragma omp taskwait

			result.add(result_temp);
		}
		else
		{
			// Terminates dimension recursion
			DataVector temp(alpha.getSize());

			#pragma omp task shared(alpha, result)
			up(alpha, result, this->algoDims[dim]);

			#pragma omp task shared(alpha, temp)
			down(alpha, temp, this->algoDims[dim]);

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
			upOpDim(alpha, temp, this->algoDims[dim]);
			updown_parallel(temp, result, dim-1, op_dim);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, op_dim);
			downOpDim(temp_two, result_temp, this->algoDims[dim]);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		upOpDim(alpha, result, this->algoDims[dim]);

		#pragma omp task shared(alpha, temp)
		downOpDim(alpha, temp, this->algoDims[dim]);

		#pragma omp taskwait

		result.add(temp);
	}
}
#endif

}
