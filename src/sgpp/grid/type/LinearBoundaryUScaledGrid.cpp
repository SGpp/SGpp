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

#include "grid/Grid.hpp"
#include "grid/type/LinearBoundaryUScaledGrid.hpp"

#include "grid/generation/BoundaryUScaledGridGenerator.hpp"

// Include all operations on the linear boundary grid
#include "basis/linearboundaryUScaled/operation/classification/OperationBLinearBoundaryUScaled.hpp"
#include "basis/linearboundaryUScaled/operation/common/OperationEvalLinearBoundaryUScaled.hpp"
#include "basis/linearboundaryUScaled/operation/common/OperationHierarchisationLinearBoundaryUScaled.hpp"
#include "basis/linearboundaryUScaled/operation/classification/OperationLaplaceLinearBoundaryUScaled.hpp"

// @todo (heinecke) removed this when done
#include "basis/linearboundaryUScaled/operation/common/OperationUpDownTestLinearBoundaryUScaled.hpp"

#include "basis/linearboundaryUScaled/operation/finance/OperationDeltaLinearTrapezoidBoundary.hpp"
#include "basis/linearboundaryUScaled/operation/finance/OperationGammaPartOneLinearTrapezoidBoundary.hpp"
#include "basis/linearboundaryUScaled/operation/finance/OperationGammaPartTwoLinearTrapezoidBoundary.hpp"
#include "basis/linearboundaryUScaled/operation/finance/OperationGammaPartThreeLinearTrapezoidBoundary.hpp"
#include "basis/linearboundaryUScaled/operation/finance/OperationRiskfreeRateLinearTrapezoidBoundary.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

LinearBoundaryUScaledGrid::LinearBoundaryUScaledGrid(std::istream& istr) : Grid(istr)
{

}

LinearBoundaryUScaledGrid::LinearBoundaryUScaledGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearBoundaryUScaledGrid::~LinearBoundaryUScaledGrid()
{
}

const char* LinearBoundaryUScaledGrid::getType()
{
	return "linearBoundaryUScaled";
}

Grid* LinearBoundaryUScaledGrid::unserialize(std::istream& istr)
{
	return new LinearBoundaryUScaledGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearBoundaryUScaledGrid::createGridGenerator()
{
	return new BoundaryUScaledGridGenerator(this->storage);
}

OperationB* LinearBoundaryUScaledGrid::createOperationB()
{
	return new OperationBLinearBoundaryUScaled(this->storage);
}

OperationMatrix* LinearBoundaryUScaledGrid::createOperationLaplace()
{
	return new OperationLaplaceLinearBoundaryUScaled(this->storage);
}

OperationEval* LinearBoundaryUScaledGrid::createOperationEval()
{
	return new OperationEvalLinearBoundaryUScaled(this->storage);
}

OperationHierarchisation* LinearBoundaryUScaledGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearBoundaryUScaled(this->storage);
}

// @todo (heinecke) removed this when done
OperationMatrix* LinearBoundaryUScaledGrid::createOperationUpDownTest()
{
	return new OperationUpDownTestLinearBoundaryUScaled(this->storage);
}

// finance operations
/////////////////////
OperationMatrix* LinearBoundaryUScaledGrid::createOperationDelta()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearBoundaryUScaledGrid::createOperationGammaPartOne()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearBoundaryUScaledGrid::createOperationGammaPartTwo()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearBoundaryUScaledGrid::createOperationGammaPartThree()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearBoundaryUScaledGrid::createOperationRiskfreeRate()
{
	throw factory_exception("Unsupported operation");
}

}
