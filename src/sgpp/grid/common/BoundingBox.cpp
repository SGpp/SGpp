/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/common/BoundingBox.hpp"

namespace sg
{

BoundingBox::BoundingBox(size_t dim)
{
	nDim = dim;
	dimensionBoundaries = new DimensionBoundary[nDim];
	for (size_t i = 0; i < nDim; i++)
	{
		dimensionBoundaries[i].leftBoundary = 0.0;
		dimensionBoundaries[i].rightBoundary = 1.0;
		dimensionBoundaries[i].bDirichletLeft = false;
		dimensionBoundaries[i].bDirichletRight = false;
	}
	bTrivialCube = true;
}

BoundingBox::BoundingBox(size_t dim, DimensionBoundary* boundaries)
{
	bTrivialCube = true;
	nDim = dim;
	dimensionBoundaries = new DimensionBoundary[nDim];
	for (size_t i = 0; i < nDim; i++)
	{
		dimensionBoundaries[i] = boundaries[i];
		if (dimensionBoundaries[i].leftBoundary != 0.0 || dimensionBoundaries[i].rightBoundary != 1.0)
		{
			bTrivialCube = false;
		}
	}
}

BoundingBox::BoundingBox(BoundingBox& copyBoundingBox)
{
	bTrivialCube = true;
	nDim = copyBoundingBox.getDimensions();
	dimensionBoundaries = new DimensionBoundary[nDim];
	for (size_t i = 0; i < nDim; i++)
	{
		dimensionBoundaries[i] = copyBoundingBox.getBoundary(i);
		if (dimensionBoundaries[i].leftBoundary != 0.0 || dimensionBoundaries[i].rightBoundary != 1.0)
		{
			bTrivialCube = false;
		}
	}
}

BoundingBox::~BoundingBox()
{
	delete[] dimensionBoundaries;
}

void BoundingBox::setBoundary(size_t dimension, DimensionBoundary& newBoundaries)
{
	dimensionBoundaries[dimension] = newBoundaries;

	if (dimensionBoundaries[dimension].leftBoundary != 0.0 || dimensionBoundaries[dimension].rightBoundary != 1.0)
	{
		bTrivialCube = false;
	}

}

DimensionBoundary BoundingBox::getBoundary(size_t dimension)
{
	return dimensionBoundaries[dimension];
}

size_t BoundingBox::getDimensions()
{
	return nDim;
}

double BoundingBox::getIntervalWidth(size_t dimension)
{
	return dimensionBoundaries[dimension].rightBoundary - dimensionBoundaries[dimension].leftBoundary;
}

double BoundingBox::getIntervalOffset(size_t dimension)
{
	return dimensionBoundaries[dimension].leftBoundary;
}

bool BoundingBox::isTrivialCube()
{
	return bTrivialCube;
}

bool BoundingBox::hasDirichletBoundaryLeft(size_t dimension)
{
	return dimensionBoundaries[dimension].bDirichletLeft;
}

bool BoundingBox::hasDirichletBoundaryRight(size_t dimension)
{
	return dimensionBoundaries[dimension].bDirichletRight;
}

}
