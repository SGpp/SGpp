/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
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
#include "grid/type/LinearGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the linear grid
#include "basis/linear/operation/datadriven/OperationBLinear.hpp"
#include "basis/linear/operation/datadriven/OperationTestLinear.hpp"
#include "basis/linear/operation/common/OperationEvalLinear.hpp"
#include "basis/linear/operation/common/OperationHierarchisationLinear.hpp"
#include "basis/linear/operation/datadriven/OperationLaplaceLinear.hpp"
#include "basis/linear/operation/pde/OperationLTwoDotProductLinear.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

LinearGrid::LinearGrid(std::istream& istr) : Grid(istr)
{

}

LinearGrid::LinearGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearGrid::~LinearGrid()
{
}

const char* LinearGrid::getType()
{
	return "linear";
}

Grid* LinearGrid::unserialize(std::istream& istr)
{
	return new LinearGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationB* LinearGrid::createOperationB()
{
	return new OperationBLinear(this->storage);
}

OperationMatrix* LinearGrid::createOperationLaplace()
{
	return new OperationLaplaceLinear(this->storage);
}

OperationEval* LinearGrid::createOperationEval()
{
	return new OperationEvalLinear(this->storage);
}

OperationTest* LinearGrid::createOperationTest()
{
	return new OperationTestLinear(this->storage);
}

OperationHierarchisation* LinearGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinear(this->storage);
}

OperationMatrix* LinearGrid::createOperationLTwoDotProduct()
{
	return new OperationLTwoDotProductLinear(this->storage);
}

// @todo (heinecke) removed this when done
OperationMatrix* LinearGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* LinearGrid::createOperationDelta(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearGrid::createOperationGamma(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

}
