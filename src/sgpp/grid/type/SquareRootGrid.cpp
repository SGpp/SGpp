/*
 * SquareRootGrid.cpp
 *
 *  Created on: Aug 4, 2010
 *      Author: Aliz Nagy
 */

/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/SquareRootGrid.hpp"

#include "grid/generation/SquareRootGridGenerator.hpp"
#include "basis/linear/boundary/operation/common/OperationEvalLinearBoundary.hpp"
#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"
#include "sgpp.hpp"

#include <iostream>

namespace sg
{

SquareRootGrid::SquareRootGrid(std::istream& istr) : Grid(istr)
{

}

SquareRootGrid::SquareRootGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

SquareRootGrid::SquareRootGrid(BoundingBox& BB)
{
	this->storage = new GridStorage(BB);
}

SquareRootGrid::~SquareRootGrid()
{
}

const char* SquareRootGrid::getType()
{
	return "squareRoot";
}
Grid* SquareRootGrid::unserialize(std::istream& istr)
{
	return new SquareRootGrid(istr);
}
/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* SquareRootGrid::createGridGenerator()
{
	return new SquareRootGridGenerator(this->storage);
}
OperationHierarchisation* SquareRootGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearBoundary(this->storage);
}
OperationEval* SquareRootGrid::createOperationEval()
{
	return new OperationEvalLinearBoundary(this->storage);
}

OperationConvert* SquareRootGrid::createOperationConvert()
{
	throw factory_exception("Unsupported operation");
}

}

