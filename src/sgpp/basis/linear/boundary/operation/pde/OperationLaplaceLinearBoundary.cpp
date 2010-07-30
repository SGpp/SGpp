/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "algorithm/common/sweep.hpp"

#include "grid/common/BoundingBox.hpp"

namespace sg
{

OperationLaplaceLinearBoundary::OperationLaplaceLinearBoundary(GridStorage* storage) : UpDownOneOpDim(storage)
{
}

OperationLaplaceLinearBoundary::~OperationLaplaceLinearBoundary()
{
}

void OperationLaplaceLinearBoundary::up(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiUpBBLinearBoundary func(this->storage);
	sweep<detail::PhiPhiUpBBLinearBoundary> s(func, this->storage);
	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::down(DataVector& alpha, DataVector& result, size_t dim)
{
	detail::PhiPhiDownBBLinearBoundary func(this->storage);
	sweep<detail::PhiPhiDownBBLinearBoundary> s(func, this->storage);
	s.sweep1D_Boundary(alpha, result, dim);
}

void OperationLaplaceLinearBoundary::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// Bounding Box handling
	BoundingBox* boundingBox = this->storage->getBoundingBox();
	double q = boundingBox->getIntervalWidth(dim);
	double Qqout = 1.0/(q*q);

	// init the coefficients of the ansatz functions with boundary
	result.setAll(0.0);

	if (q != 1.0)
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			if (level == 0)
			{
				if (index == 0)
				{
					if (!boundingBox->hasDirichletBoundaryLeft(dim))
					{
						//only affects the diagonal of the stiffness matrix
						result[i] += Qqout*alpha[i];

						// down
						if (index == 0)
						{
							GridIndex index_one = (*storage)[i];
							index_one.set(dim, 0, 1);
							if (!boundingBox->hasDirichletBoundaryRight(dim))
							{
								result[(*storage)[&index_one]] += ((-1.0*Qqout) * alpha[i]);
							}
						}
					}
				}
				if (index == 1)
				{
					if (!boundingBox->hasDirichletBoundaryRight(dim))
					{
						//only affects the diagonal of the stiffness matrix
						result[i] += Qqout*alpha[i];
					}
				}
			}
			//only affects the diagonal of the stiffness matrix
			else
			{
				result[i] = alpha[i]*(Qqout*pow(2.0, static_cast<int>(level+1)));
			}
		}
	}
	else
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			if (level == 0)
			{
				if (index == 0)
				{
					if (!boundingBox->hasDirichletBoundaryLeft(dim))
					{
						//only affects the diagonal of the stiffness matrix
						result[i] += alpha[i];

						// down
						if (index == 0)
						{
							GridIndex index_one = (*storage)[i];
							index_one.set(dim, 0, 1);
							if (!boundingBox->hasDirichletBoundaryRight(dim))
							{
								result[(*storage)[&index_one]] += ((-1.0) * alpha[i]);
							}
						}
					}
				}
				if (index == 1)
				{
					if (!boundingBox->hasDirichletBoundaryRight(dim))
					{
						//only affects the diagonal of the stiffness matrix
						result[i] += alpha[i];
					}
				}
			}
			//only affects the diagonal of the stiffness matrix
			else
			{
				result[i] = alpha[i]*pow(2.0, static_cast<int>(level+1));
			}
		}
	}
}

void OperationLaplaceLinearBoundary::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	// Bounding Box handling
	BoundingBox* boundingBox = this->storage->getBoundingBox();
	double q = boundingBox->getIntervalWidth(dim);
	double Qqout = 1.0/(q*q);

	// init the coefficients of the ansatz functions with boundary
	result.setAll(0.0);

	if (q != 1.0)
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			if (level == 0)
			{
				// up
				if (index == 1)
				{
					GridIndex index_zero = (*storage)[i];
					index_zero.set(dim, 0, 0);
					if (!boundingBox->hasDirichletBoundaryLeft(dim))
					{
						result[(*storage)[&index_zero]] += ((-1.0*Qqout) * alpha[i]);
					}
				}
			}
		}
	}
	else
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			if (level == 0)
			{
				// up
				if (index == 1)
				{
					GridIndex index_zero = (*storage)[i];
					index_zero.set(dim, 0, 0);
					if (!boundingBox->hasDirichletBoundaryLeft(dim))
					{
						result[(*storage)[&index_zero]] += ((-1.0) * alpha[i]);
					}
				}
			}
		}
	}
}

}
