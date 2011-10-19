/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "grid/Grid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"

#include "grid/generation/BoundaryGridGenerator.hpp"

#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"

#include "exception/factory_exception.hpp"


#include <iostream>

namespace sg
{
namespace base
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
	return new BoundaryGridGenerator(this->storage);
}

}
}
