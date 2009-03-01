/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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
#include "grid/type/LinearBoundaryOScaledGrid.hpp"

#include "grid/generation/BoundaryOScaledGridGenerator.hpp"

// Include all operations on the linear boundary grid
#include "basis/linearboundaryOScaled/operation/OperationBLinearBoundaryOScaled.hpp"
#include "basis/linearboundaryOScaled/operation/OperationEvalLinearBoundaryOScaled.hpp"
#include "basis/linearboundaryOScaled/operation/OperationHierarchisationLinearBoundaryOScaled.hpp"
#include "basis/linearboundaryOScaled/operation/OperationLaplaceLinearBoundaryOScaled.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

LinearBoundaryOScaledGrid::LinearBoundaryOScaledGrid(std::istream& istr) : Grid(istr)
{

}

LinearBoundaryOScaledGrid::LinearBoundaryOScaledGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearBoundaryOScaledGrid::~LinearBoundaryOScaledGrid()
{
}

const char* LinearBoundaryOScaledGrid::getType()
{
	return "linearBoundaryOScaled";
}

Grid* LinearBoundaryOScaledGrid::unserialize(std::istream& istr)
{
	return new LinearBoundaryOScaledGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearBoundaryOScaledGrid::createGridGenerator()
{
	return new BoundaryOScaledGridGenerator(this->storage);
}

OperationB* LinearBoundaryOScaledGrid::createOperationB()
{
	return new OperationBLinearBoundaryOScaled(this->storage);
}

OperationMatrix* LinearBoundaryOScaledGrid::createOperationLaplace()
{
	return new OperationLaplaceLinearBoundaryOScaled(this->storage);
}

OperationEval* LinearBoundaryOScaledGrid::createOperationEval()
{
	return new OperationEvalLinearBoundaryOScaled(this->storage);
}

OperationHierarchisation* LinearBoundaryOScaledGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearBoundaryOScaled(this->storage);
}

}
