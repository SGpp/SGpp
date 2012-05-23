/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/algorithm/UpDownFourOpDims.hpp"

namespace sg
{
namespace pde
{

UpDownFourOpDims::UpDownFourOpDims(sg::base::GridStorage* storage, double**** coef) : storage(storage), coefs(coef), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size())
{
}

UpDownFourOpDims::UpDownFourOpDims(sg::base::GridStorage* storage) : storage(storage), coefs(NULL), algoDims(storage->getAlgorithmicDimensions()), numAlgoDims_(storage->getAlgorithmicDimensions().size())
{
}

UpDownFourOpDims::~UpDownFourOpDims()
{
}

void UpDownFourOpDims::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
{
	result.setAll(0.0);

#pragma omp parallel
	{
#pragma omp single nowait
		{
			for(size_t i = 0; i < this->numAlgoDims_; i++)
			{
				for(size_t j = 0; j < this->numAlgoDims_; j++)
				{
					for(size_t k = 0; k < this->numAlgoDims_; k++)
					{
						for(size_t l = 0; l < this->numAlgoDims_; l++)
						{

#pragma omp task firstprivate(i, j) shared(alpha, result)
							{
								sg::base::DataVector beta(result.getSize());

								if (this->coefs != NULL)
								{
									if (this->coefs[i][j][k][l] != 0.0)
									{
										this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j,k,l);

#pragma omp critical
										{
											result.axpy(this->coefs[i][j][k][l],beta);
										}
									}
								}
								else
								{
									this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j,k,l);

#pragma omp critical
									{
										result.add(beta);
									}
								}
							}

						}
					}
				}
			}

#pragma omp taskwait
		}
	}
}

void UpDownFourOpDims::specialOpX(sg::base::DataVector& alpha, sg::base::DataVector& result, void (sg::pde::UpDownFourOpDims::*pt2UpFunc)(sg::base::DataVector&, sg::base::DataVector&, size_t), void (sg::pde::UpDownFourOpDims::*pt2DownFunc)(sg::base::DataVector&, sg::base::DataVector&, size_t), size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		sg::base::DataVector temp(alpha.getSize());
		sg::base::DataVector result_temp(alpha.getSize());
		sg::base::DataVector temp_two(alpha.getSize());

#pragma omp task if(this->numAlgoDims_ - dim <= this->maxParallelDims_) shared(alpha, temp, result)
		{
			(this->*pt2UpFunc)(alpha, temp, this->algoDims[dim]);
			updown(temp, result, dim-1, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
		}

		// Same from the other direction:
#pragma omp task if(this->numAlgoDims_ - dim <= this->maxParallelDims_) shared(alpha, temp_two, result_temp)
		{
			updown(alpha, temp_two, dim-1, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
			(this->*pt2DownFunc)(temp_two, result_temp, this->algoDims[dim]);
		}

#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		sg::base::DataVector temp(alpha.getSize());

#pragma omp task if(this->numAlgoDims_ - dim <= this->maxParallelDims_) shared(alpha, result)
		(this->*pt2UpFunc)(alpha, result, this->algoDims[dim]);

#pragma omp task if(this->numAlgoDims_ - dim <= this->maxParallelDims_) shared(alpha, temp)
		(this->*pt2DownFunc)(alpha, temp, this->algoDims[dim]);

#pragma omp taskwait

		result.add(temp);
	}
}


void UpDownFourOpDims::specialOpOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four)
{
	specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimOne, &sg::pde::UpDownFourOpDims::downOpDimOne, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four)
{
	specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimTwo, &sg::pde::UpDownFourOpDims::downOpDimTwo, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);

}

void UpDownFourOpDims::specialOpThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four)
{
	specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimThree, &sg::pde::UpDownFourOpDims::downOpDimThree, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

void UpDownFourOpDims::specialOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four)
{
	specialOpX(alpha, result, &sg::pde::UpDownFourOpDims::upOpDimFour, &sg::pde::UpDownFourOpDims::downOpDimFour, dim, op_dim_one, op_dim_two, op_dim_three, op_dim_four);
}

}
}
