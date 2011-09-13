/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#include "grid/Grid.hpp"
#include "grid/type/PolyGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

#include "exception/factory_exception.hpp"


#include <iostream>

namespace sg
{
namespace base
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

size_t PolyGrid::getDegree() const
{
	return this->degree;
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

}
}
