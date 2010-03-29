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

#include "grid/Grid.hpp"
#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"

#include "grid/generation/TrapezoidBoundaryGridGenerator.hpp"

// Include all operations on the linear boundary grid
#include "basis/linear/boundary/operation/datadriven/OperationBLinearBoundary.hpp"
#include "basis/linear/boundary/operation/datadriven/OperationTestLinearBoundary.hpp"
#include "basis/linear/boundary/operation/common/OperationEvalLinearBoundary.hpp"
#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"
// @todo (heinecke) removed this when done
#include "basis/linear/boundary/operation/common/OperationUpDownTestLinearBoundary.hpp"

#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/OperationLTwoDotProductLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLinearBoundary.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(std::istream& istr) : Grid(istr)
{

}

LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(BoundingBox& BB)
{
	this->storage = new GridStorage(BB);
}

LinearTrapezoidBoundaryGrid::~LinearTrapezoidBoundaryGrid()
{
}

const char* LinearTrapezoidBoundaryGrid::getType()
{
	return "linearTrapezoidBoundary";
}

Grid* LinearTrapezoidBoundaryGrid::unserialize(std::istream& istr)
{
	return new LinearTrapezoidBoundaryGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearTrapezoidBoundaryGrid::createGridGenerator()
{
	return new TrapezoidBoundaryGridGenerator(this->storage);
}

OperationB* LinearTrapezoidBoundaryGrid::createOperationB()
{
	return new OperationBLinearBoundary(this->storage);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationLaplace()
{
	return new OperationLaplaceLinearBoundary(this->storage);
}

OperationEval* LinearTrapezoidBoundaryGrid::createOperationEval()
{
	return new OperationEvalLinearBoundary(this->storage);
}

OperationTest* LinearTrapezoidBoundaryGrid::createOperationTest()
{
	return new OperationTestLinearBoundary(this->storage);
}

OperationHierarchisation* LinearTrapezoidBoundaryGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearBoundary(this->storage);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationLTwoDotProduct()
{
	return new OperationLTwoDotProductLinearBoundary(this->storage);
}

// @todo (heinecke) removed this when done
OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationUpDownTest()
{
	return new OperationUpDownTestLinearBoundary(this->storage);
}

// finance operations
/////////////////////
OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationDelta(DataVector& coef)
{
	return new OperationDeltaLinearBoundary(this->storage, coef);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationGamma(DataVector& coef)
{
	return new OperationGammaLinearBoundary(this->storage, coef);
}

}
