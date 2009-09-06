/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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
	dimensionBoundaries[dimension].leftBoundary = newBoundaries.leftBoundary;
	dimensionBoundaries[dimension].rightBoundary = newBoundaries.rightBoundary;
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

}
