/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
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
#include "grid/type/PolyGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the poly grid
#include "basis/poly/operation/OperationBPoly.hpp"
#include "basis/poly/operation/OperationEvalPoly.hpp"
#include "basis/poly/operation/OperationHierarchisationPoly.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

PolyGrid::PolyGrid(std::istream& istr) : Grid(istr), degree(-1)
{
	istr >> degree;
}

PolyGrid::PolyGrid(size_t dim, size_t degree) : degree(degree)
{
	this->storage = new GridStorage(dim);
}

PolyGrid::~PolyGrid()
{
}

const char* PolyGrid::getType()
{
	return "poly";
}

Grid* PolyGrid::unserialize(std::istream& istr)
{
	return new PolyGrid(istr);
}

void PolyGrid::serialize(std::ostream& ostr)
{
	this->Grid::serialize(ostr);
	ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* PolyGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationB* PolyGrid::createOperationB()
{
	return new OperationBPoly(this->storage, this->degree);
}

OperationMatrix* PolyGrid::createOperationLaplace()
{
	throw factory_exception("Unsupported operation");
}

OperationEval* PolyGrid::createOperationEval()
{
	return new OperationEvalPoly(this->storage, this->degree);
}

OperationHierarchisation* PolyGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationPoly(this->storage, this->degree);
}

}
