/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2008  Joerg Blank (blankj@in.tum.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "GridFactory.hpp"
#include "GridFactoryImpl.hpp"

#include "Operations.hpp"
#include "StandardOperations.hpp"
#include "linear/LinearOperations.hpp"
#include "modlinear/ModLinearOperations.hpp"
#include "poly/PolyOperations.hpp"
#include "modpoly/ModPolyOperations.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

// ***** LinearGrid ***** //

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

OperationHierarchisation* LinearGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinear(this->storage);
}

// ***** ModLinearGrid *****

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

OperationHierarchisation* ModLinearGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationModLinear(this->storage);
}

// ***** PolyGrid *****

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

// ***** ModPolyGrid ***** //

ModPolyGrid::ModPolyGrid(std::istream& istr) : Grid(istr), degree(-1)
{
	istr >> degree;
}


ModPolyGrid::ModPolyGrid(size_t dim, size_t degree) : degree(degree)
{
	this->storage = new GridStorage(dim);
}

ModPolyGrid::~ModPolyGrid()
{
}

const char* ModPolyGrid::getType()
{
	return "modpoly";
}

Grid* ModPolyGrid::unserialize(std::istream& istr)
{
	return new ModPolyGrid(istr);
}

void ModPolyGrid::serialize(std::ostream& ostr)
{
	this->Grid::serialize(ostr);
	ostr << degree << std::endl;
}


/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* ModPolyGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationB* ModPolyGrid::createOperationB()
{
	return new OperationBModPoly(this->storage, this->degree);
}

OperationMatrix* ModPolyGrid::createOperationLaplace()
{
	throw factory_exception("Unsupported operation");
}

OperationEval* ModPolyGrid::createOperationEval()
{
	return new OperationEvalModPoly(this->storage, this->degree);
}

OperationHierarchisation* ModPolyGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationModPoly(this->storage, this->degree);
}

}
