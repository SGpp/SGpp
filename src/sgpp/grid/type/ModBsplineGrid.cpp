/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2010 Dirk Pflueger (pflueged@in.tum.de)                     */
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
#include "grid/type/ModBsplineGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the mod bspline grid
#include "basis/modbspline/operation/datadriven/OperationBModBspline.hpp"
#include "basis/modbspline/operation/datadriven/OperationTestModBspline.hpp"
#include "basis/modbspline/operation/common/OperationEvalModBspline.hpp"
#include "basis/modbspline/operation/common/OperationHierarchisationModBspline.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

ModBsplineGrid::ModBsplineGrid(std::istream& istr) : Grid(istr), degree(-1)
{
    istr >> degree;
}


ModBsplineGrid::ModBsplineGrid(size_t dim, size_t degree) : degree(degree)
{
    this->storage = new GridStorage(dim);
}

ModBsplineGrid::~ModBsplineGrid()
{
}

const char* ModBsplineGrid::getType()
{
	return "modBspline";
}

Grid* ModBsplineGrid::unserialize(std::istream& istr)
{
	return new ModBsplineGrid(istr);
}

void ModBsplineGrid::serialize(std::ostream& ostr)
{
    this->Grid::serialize(ostr);
    ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* ModBsplineGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationB* ModBsplineGrid::createOperationB()
{
	return new OperationBModBspline(this->storage, this->degree);
}

OperationMatrix* ModBsplineGrid::createOperationLaplace()
{
	throw factory_exception("Unsupported operation");
}

OperationEval* ModBsplineGrid::createOperationEval()
{
	return new OperationEvalModBspline(this->storage, this->degree);
}

OperationTest* ModBsplineGrid::createOperationTest()
{
    return new OperationTestModBspline(this->storage, this->degree);
}

OperationHierarchisation* ModBsplineGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationModBspline(this->storage, this->degree);
}

OperationMatrix* ModBsplineGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModBsplineGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* ModBsplineGrid::createOperationDelta(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModBsplineGrid::createOperationGamma(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

}
