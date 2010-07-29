/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/ModLinearGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the mod linear grid
#include "basis/modlinear/operation/datadriven/OperationBModLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationTestModLinear.hpp"
#include "basis/modlinear/operation/common/OperationEvalModLinear.hpp"
#include "basis/modlinear/operation/common/OperationHierarchisationModLinear.hpp"

#include "exception/factory_exception.hpp"

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


}
