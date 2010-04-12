/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#include "grid/common/DirichletGridConverter.hpp"
#include "grid/Grid.hpp"
#include "grid/type/LinearGrid.hpp"

#include "exception/generation_exception.hpp"

#include <cstring>

namespace sg
{

DirichletGridConverter::DirichletGridConverter() : numTotalGridPoints(0), numInnerGridPoints(0), conCoefArray(NULL), bFirstTime(false)
{
}

DirichletGridConverter::~DirichletGridConverter()
{
	if (this->numInnerGridPoints > 0)
	{
		delete[] this->conCoefArray;
	}
}

void DirichletGridConverter::buildInnerGridWithCoefs(Grid& BoundaryGrid, DataVector& BoundaryCoefs, Grid* InnerGrid, DataVector* InnerCoefs)
{
	if (this->bFirstTime == true)
	{
		if (strcmp(BoundaryGrid.getType(), "linearBoundary") == 0 || strcmp(BoundaryGrid.getType(), "linearTrapezoidBoundary") == 0)
		{
			GridStorage* myGridStorage = BoundaryGrid.getStorage();

			// determine the number of grid points for both grids
			this->numTotalGridPoints = myGridStorage->size();
			this->numInnerGridPoints = myGridStorage->getNumInnerPoints();

			// allocate the translation array for the coefficients
			this->conCoefArray = new size_t[this->numInnerGridPoints];

			// create new inner Grid, with one grid point
			InnerGrid = new LinearGrid(*BoundaryGrid.getBoundingBox());

			// create new DataVector for storing the inner grid's coefficients
			InnerCoefs = new DataVector(this->numInnerGridPoints);

			// Iterate through all grid points and filter inner points
			size_t numInner = 0;
			for (size_t i = 0; i < this->numTotalGridPoints; i++)
			{
				GridIndex* curPoint = (*myGridStorage)[i];
				if (curPoint->isInnerPoint() == true)
				{
					// handle coefficients
					this->conCoefArray[numInner] = i;
					InnerCoefs[numInner] = BoundaryCoefs[i];
					numInner++;
					// insert point into inner grid
					InnerGrid->getStorage()->insert(*curPoint);
				}
			}
			this->bFirstTime = true;
		}
		else
		{
			throw generation_exception("DirichletGridConverter : buildInnerGridWithCoefs : Boundary Grid is from an unsupported grid type!");
		}
	}
	else
	{
		throw generation_exception("DirichletGridConverter : buildInnerGridWithCoefs : This method can only be called once for one instance!");
	}
}

void DirichletGridConverter::rebuildInnerGridWithCoefs(Grid& BoundaryGrid, DataVector& BoundaryCoefs, Grid* InnerGrid, DataVector* InnerCoefs)
{
	if (this->bFirstTime == false)
	{
		if (strcmp(BoundaryGrid.getType(), "linearBoundary") == 0 || strcmp(BoundaryGrid.getType(), "linearTrapezoidBoundary") == 0)
		{
			GridStorage* myGridStorage = BoundaryGrid.getStorage();

			// determine the number of grid points for both grids
			this->numTotalGridPoints = myGridStorage->size();
			this->numInnerGridPoints = myGridStorage->getNumInnerPoints();

			// allocate the translation array for the coefficients
			delete[] this->conCoefArray;
			this->conCoefArray = new size_t[this->numInnerGridPoints];

			// create new inner Grid, with one grid point
			InnerGrid->getStorage()->emptyStorage();

			// create new DataVector for storing the inner grid's coefficients
			delete[] InnerCoefs;
			InnerCoefs = new DataVector(this->numInnerGridPoints);

			// Iterate through all grid points and filter inner points
			size_t numInner = 0;
			for (size_t i = 0; i < this->numTotalGridPoints; i++)
			{
				GridIndex* curPoint = (*myGridStorage)[i];
				if (curPoint->isInnerPoint() == true)
				{
					// handle coefficients
					this->conCoefArray[numInner] = i;
					InnerCoefs[numInner] = BoundaryCoefs[i];
					numInner++;
					// insert point into inner grid
					InnerGrid->getStorage()->insert(*curPoint);
				}
			}
		}
		else
		{
			throw generation_exception("DirichletGridConverter : rebuildInnerGridWithCoefs : Boundary Grid is from an unsupported grid type!");
		}
	}
	else
	{
		throw generation_exception("DirichletGridConverter : rebuildInnerGridWithCoefs : This method can only be called after initial build!");
	}
}

void DirichletGridConverter::calcInnerCoefs(DataVector& BoundaryCoefs, DataVector& InnerCoefs)
{
	if (this->numInnerGridPoints == 0)
	{
		for (size_t i = 0; i < this->numInnerGridPoints; i++)
		{
			InnerCoefs[i] = BoundaryCoefs[this->conCoefArray[i]];
		}
	}
	else
	{
		throw generation_exception("DirichletGridConverter : calcInnerCoefs : Inner grid  hasn't be created before!");
	}
}

void DirichletGridConverter::updateBoundaryCoefs(DataVector& BoundaryCoefs, DataVector& InnerCoefs)
{
	if (this->numInnerGridPoints == 0)
	{
		for (size_t i = 0; i < this->numInnerGridPoints; i++)
		{
			BoundaryCoefs[this->conCoefArray[i]] = InnerCoefs[i];
		}
	}
	else
	{
		throw generation_exception("DirichletGridConverter : calcInnerCoefs : Inner grid  hasn't be created before!");
	}
}

}
