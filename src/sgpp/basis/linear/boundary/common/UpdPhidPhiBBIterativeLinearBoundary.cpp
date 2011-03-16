/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de

#include "basis/linear/boundary/common/UpdPhidPhiBBIterativeLinearBoundary.hpp"
#include "grid/common/BoundingBox.hpp"
using namespace sg::base;

namespace sg
{

UpdPhidPhiBBIterativeLinearBoundary::UpdPhidPhiBBIterativeLinearBoundary(GridStorage* storage) : storage(storage)
{
}

UpdPhidPhiBBIterativeLinearBoundary::~UpdPhidPhiBBIterativeLinearBoundary()
{
}

void UpdPhidPhiBBIterativeLinearBoundary::operator()(DataVector& alpha, DataVector& result, size_t dim)
{
	// Bounding Box handling
	BoundingBox* boundingBox = this->storage->getBoundingBox();
	double q = boundingBox->getIntervalWidth(dim);
	double Qqout = 1.0/q;

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
