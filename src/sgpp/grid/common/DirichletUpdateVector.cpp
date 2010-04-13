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

void DirichletUpdateVector::setInnerPointsToZero(DataVector& updateVector)
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
			if (level > 0)
			{
				if (myBounds.bDirichletRight == true && myBounds.bDirichletLeft == true)
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
