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
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == false)
		{
			updateVector.set(i, sourceVector.get(i));
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
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == false)
		{
			updateVector.set(i, 0.0);
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
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == true)
		{
			updateVector.set(i, 0.0);
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
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == false)
		{
			updateVector.set(i, updateVector.get(i)*value);
		}
	}
}

}
