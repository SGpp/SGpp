/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/UpDownTwoOpDims.hpp"

namespace sg
{

UpDownTwoOpDims::UpDownTwoOpDims(GridStorage* storage, DataVector& coef)
{
	this->storage = storage;
	this->coefs = &coef;
}

UpDownTwoOpDims::UpDownTwoOpDims(GridStorage* storage)
{
	this->storage = storage;
	this->coefs = NULL;
}


UpDownTwoOpDims::~UpDownTwoOpDims()
{
}

void UpDownTwoOpDims::mult(DataVector& alpha, DataVector& result)
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
#ifndef AIX_XLC
					#pragma omp single nowait
#endif
#ifdef AIX_XLC
					#pragma omp single
#endif
					{
						if (this->coefs != NULL)
						{
							if (this->coefs->get((storage->dim()*i)+j) != 0.0)
							{
								this->updown_parallel(alpha, beta, storage->dim() - 1, i, j);
							}
						}
						else
						{
							this->updown_parallel(alpha, beta, storage->dim() - 1, i, j);
						}
					}
				}
			}
			// Calculate the "diagonal" of the operation
			if (j <= i)
			{
				if (this->coefs != NULL)
				{
					result.axpy(this->coefs->get((storage->dim()*i)+j),beta);
				}
				else
				{
					result.add(beta);
				}
			}
		}
	}
#endif
#ifndef USEOMPTHREE
	#pragma omp parallel
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
					if (this->coefs != NULL)
					{
						if (this->coefs->get((storage->dim()*i)+j) != 0.0)
						{
							this->updown(alpha, beta, storage->dim() - 1, i, j);
						}
					}
					else
					{
						this->updown(alpha, beta, storage->dim() - 1, i, j);
					}
				}

				// Calculate the "diagonal" of the operation
				if (j <= i)
				{
					#pragma omp critical
					{
						if (this->coefs != NULL)
						{
							result.axpy(this->coefs->get((storage->dim()*i)+j),beta);
						}
						else
						{
							result.add(beta);
						}
					}
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
			if (j <= i)
			{
				if (this->coefs != NULL)
				{
					result.axpy(this->coefs->get((storage->dim()*i)+j),beta);
				}
				else
				{
					result.add(beta);
				}
			}
		}
	}
#endif
}

#ifndef USEOMPTHREE
void UpDownTwoOpDims::updown(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two)
{
	if((dim == op_dim_one) && (dim == op_dim_two))
	{
		specialOpOneAndOpTwo(alpha, result, dim, op_dim_one, op_dim_two);
	}
	else if ((dim == op_dim_one || dim == op_dim_two) && (op_dim_one != op_dim_two))
	{
		if (dim == op_dim_one)
		{
			specialOpOne(alpha, result, dim, op_dim_one, op_dim_two);
		}
		if (dim == op_dim_two)
		{
			specialOpTwo(alpha, result, dim, op_dim_one, op_dim_two);
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
			updown(temp, result, dim-1, op_dim_one, op_dim_two);


			// Same from the other direction:
			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1, op_dim_one, op_dim_two);
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

void UpDownTwoOpDims::specialOpOne(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upOpDimOne(alpha, temp, dim);
		updown(temp, result, dim-1, op_dim_one, op_dim_two);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, op_dim_one, op_dim_two);
		downOpDimOne(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		upOpDimOne(alpha, result, dim);

		DataVector temp(alpha.getSize());
		downOpDimOne(alpha, temp, dim);

		result.add(temp);
	}
}

void UpDownTwoOpDims::specialOpTwo(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upOpDimTwo(alpha, temp, dim);
		updown(temp, result, dim-1, op_dim_one, op_dim_two);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, op_dim_one, op_dim_two);
		downOpDimTwo(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		upOpDimTwo(alpha, result, dim);

		DataVector temp(alpha.getSize());
		downOpDimTwo(alpha, temp, dim);

		result.add(temp);
	}
}

void UpDownTwoOpDims::specialOpOneAndOpTwo(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		upOpDimOneAndOpDimTwo(alpha, temp, dim);
		updown(temp, result, dim-1, op_dim_one, op_dim_two);


		// Same from the other direction:
		DataVector result_temp(alpha.getSize());
		updown(alpha, temp, dim-1, op_dim_one, op_dim_two);
		downOpDimOneAndOpDimTwo(temp, result_temp, dim);

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		upOpDimOneAndOpDimTwo(alpha, result, dim);

		DataVector temp(alpha.getSize());
		downOpDimOneAndOpDimTwo(alpha, temp, dim);

		result.add(temp);
	}
}
#endif

#ifdef USEOMPTHREE
void UpDownTwoOpDims::updown_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two)
{
	if((dim == op_dim_one) && (dim == op_dim_two))
	{
		specialOpOneAndOpTwo_parallel(alpha, result, dim, op_dim_one, op_dim_two);
	}
	else if ((dim == op_dim_one || dim == op_dim_two) && (op_dim_one != op_dim_two))
	{
		if (dim == op_dim_one)
		{
			specialOpOne_parallel(alpha, result, dim, op_dim_one, op_dim_two);
		}
		if (dim == op_dim_two)
		{
			specialOpTwo_parallel(alpha, result, dim, op_dim_one, op_dim_two);
		}
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
				updown_parallel(temp, result, dim-1, op_dim_one, op_dim_two);
			}

			// Same from the other direction:
			#pragma omp task shared(alpha, temp_two, result_temp)
			{
				updown_parallel(alpha, temp_two, dim-1, op_dim_one, op_dim_two);
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

void UpDownTwoOpDims::specialOpOne_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two)
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
			upOpDimOne(alpha, temp, dim);
			updown_parallel(temp, result, dim-1, op_dim_one, op_dim_two);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, op_dim_one, op_dim_two);
			downOpDimOne(temp_two, result_temp, dim);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		upOpDimOne(alpha, result, dim);

		#pragma omp task shared(alpha, temp)
		downOpDimOne(alpha, temp, dim);

		#pragma omp taskwait

		result.add(temp);
	}
}

void UpDownTwoOpDims::specialOpTwo_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two)
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
			upOpDimTwo(alpha, temp, dim);
			updown_parallel(temp, result, dim-1, op_dim_one, op_dim_two);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, op_dim_one, op_dim_two);
			downOpDimTwo(temp_two, result_temp, dim);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		upOpDimTwo(alpha, result, dim);

		#pragma omp task shared(alpha, temp)
		downOpDimTwo(alpha, temp, dim);

		#pragma omp taskwait

		result.add(temp);
	}
}

void UpDownTwoOpDims::specialOpOneAndOpTwo_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two)
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
			upOpDimOneAndOpDimTwo(alpha, temp, dim);
			updown_parallel(temp, result, dim-1, op_dim_one, op_dim_two);
		}

		// Same from the other direction:
		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown_parallel(alpha, temp_two, dim-1, op_dim_one, op_dim_two);
			downOpDimOneAndOpDimTwo(temp_two, result_temp, dim);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		upOpDimOneAndOpDimTwo(alpha, result, dim);

		#pragma omp task shared(alpha, temp)
		downOpDimOneAndOpDimTwo(alpha, temp, dim);

		#pragma omp taskwait

		result.add(temp);
	}
}
#endif

}
