
/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/LinearStretchedGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the linearstretched grid
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationTestLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/common/OperationEvalLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/common/OperationHierarchisationLinearStretched.hpp"

#include "basis/linearstretched/noboundary/operation/pde/OperationLaplaceLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/OperationLTwoDotProductLinearStretched.hpp"

#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLogLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLogLinearStretched.hpp"

#include "exception/factory_exception.hpp"


#include <iostream>

namespace sg
{
namespace base
{

LinearStretchedGrid::LinearStretchedGrid(std::istream& istr) : Grid(istr)
{

}

LinearStretchedGrid::LinearStretchedGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearStretchedGrid::LinearStretchedGrid(Stretching& BB)
{
	this->storage = new GridStorage(BB);
}

LinearStretchedGrid::~LinearStretchedGrid()
{
}

const char* LinearStretchedGrid::getType()
{
	return "linearStretched";
}

Grid* LinearStretchedGrid::unserialize(std::istream& istr)
{
	return new LinearStretchedGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearStretchedGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}
}
}
