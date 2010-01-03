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
#include "grid/type/ModLinearGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the mod linear grid
#include "basis/modlinear/operation/datadriven/OperationBModLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationTestModLinear.hpp"
#include "basis/modlinear/operation/common/OperationEvalModLinear.hpp"
#include "basis/modlinear/operation/common/OperationHierarchisationModLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationLaplaceModLinear.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

ModLinearGrid::ModLinearGrid(std::istream& istr) : Grid(istr)
{
}

ModLinearGrid::ModLinearGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

ModLinearGrid::~ModLinearGrid()
{
}

const char* ModLinearGrid::getType()
{
	return "modlinear";
}

Grid* ModLinearGrid::unserialize(std::istream& istr)
{
	return new ModLinearGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* ModLinearGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationB* ModLinearGrid::createOperationB()
{
	return new OperationBModLinear(this->storage);
}

OperationMatrix* ModLinearGrid::createOperationLaplace()
{
	return new OperationLaplaceModLinear(this->storage);
}

OperationEval* ModLinearGrid::createOperationEval()
{
	return new OperationEvalModLinear(this->storage);
}

OperationTest* ModLinearGrid::createOperationTest()
{
	return new OperationTestModLinear(this->storage);
}

OperationHierarchisation* ModLinearGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationModLinear(this->storage);
}

OperationMatrix* ModLinearGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}

// @todo (heinecke) removed this when done
OperationMatrix* ModLinearGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* ModLinearGrid::createOperationDelta(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModLinearGrid::createOperationGamma(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

}
