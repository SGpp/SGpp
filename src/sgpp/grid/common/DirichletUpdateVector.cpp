/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/common/DirichletUpdateVector.hpp"

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg
{

DirichletUpdateVector::DirichletUpdateVector(GridStorage* storage): myBoundingBox(storage->getBoundingBox()), storage(storage)
{
}

DirichletUpdateVector::~DirichletUpdateVector()
{
}

void DirichletUpdateVector::applyDirichletConditions(DataVector& updateVector, DataVector& sourceVector)
{
#ifdef USEOMP
	#pragma omp parallel for shared(updateVector) schedule(static)
#endif
	for (size_t i = 0; i < storage->size(); i++)
	{
		GridStorage::index_type::level_type level;
		GridStorage::index_type::index_type index;
		size_t dim;
		DimensionBoundary myBounds;
		for (size_t j = 0; j < storage->dim(); j++)
		{
			dim = j;
			(*storage)[i]->get(dim, level, index);
			myBounds = myBoundingBox->getBoundary(dim);
			if (level == 0)
			{
				if (index == 0 && myBounds.bDirichletLeft == true)
				{
					updateVector.set(i, sourceVector.get(i));
				}
				if (index == 1 && myBounds.bDirichletRight == true)
				{
					updateVector.set(i, sourceVector.get(i));
				}
			}
		}
	}
}

void DirichletUpdateVector::setBoundariesToZero(DataVector& updateVector)
{
#ifdef USEOMP
	#pragma omp parallel for shared(updateVector) schedule(static)
#endif
	for (size_t i = 0; i < storage->size(); i++)
	{
		GridStorage::index_type::level_type level;
		GridStorage::index_type::index_type index;
		size_t dim;
		DimensionBoundary myBounds;
		for (size_t j = 0; j < storage->dim(); j++)
		{
			dim = j;
			(*storage)[i]->get(dim, level, index);
			myBounds = myBoundingBox->getBoundary(dim);
			if (level == 0)
			{
				if (index == 0 && myBounds.bDirichletLeft == true)
				{
					updateVector.set(i, 0.0);
				}
				if (index == 1 && myBounds.bDirichletRight == true)
				{
					updateVector.set(i, 0.0);
				}
			}
		}
	}
}

void DirichletUpdateVector::multiplyBoundary(DataVector& updateVector, double value)
{
#ifdef USEOMP
	#pragma omp parallel for shared(updateVector) schedule(static)
#endif
	for (size_t i = 0; i < storage->size(); i++)
	{
		GridStorage::index_type::level_type level;
		GridStorage::index_type::index_type index;
		size_t dim;
		DimensionBoundary myBounds;
		for (size_t j = 0; j < storage->dim(); j++)
		{
			dim = j;
			(*storage)[i]->get(dim, level, index);
			myBounds = myBoundingBox->getBoundary(dim);
			if (level == 0)
			{
				if (index == 0 && myBounds.bDirichletLeft == true)
				{
					updateVector.set(i, updateVector.get(i)*value);
				}
				if (index == 1 && myBounds.bDirichletRight == true)
				{
					updateVector.set(i, updateVector.get(i)*value);
				}
			}
		}
	}
}

}
