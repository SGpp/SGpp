/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#include "grid/Grid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the linear boundary grid
#include "basis/linearboundary/operation/OperationBLinearBoundary.hpp"
#include "basis/linearboundary/operation/OperationEvalLinearBoundary.hpp"
#include "basis/linearboundary/operation/OperationHierarchisationLinearBoundary.hpp"
#include "basis/linearboundary/operation/OperationLaplaceLinearBoundary.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

LinearBoundaryGrid::LinearBoundaryGrid(std::istream& istr) : Grid(istr)
{

}

LinearBoundaryGrid::LinearBoundaryGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearBoundaryGrid::~LinearBoundaryGrid()
{
}

const char* LinearBoundaryGrid::getType()
{
	return "linearBoundary";
}

Grid* LinearBoundaryGrid::unserialize(std::istream& istr)
{
	return new LinearBoundaryGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearBoundaryGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationB* LinearBoundaryGrid::createOperationB()
{
	return NULL;
}

OperationMatrix* LinearBoundaryGrid::createOperationLaplace()
{
	return NULL;
}

OperationEval* LinearBoundaryGrid::createOperationEval()
{
	return NULL;
}

OperationHierarchisation* LinearBoundaryGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearBoundary(this->storage);
}

}
