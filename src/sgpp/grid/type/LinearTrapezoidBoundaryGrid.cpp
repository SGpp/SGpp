/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"

#include "grid/generation/TrapezoidBoundaryGridGenerator.hpp"

// Include all operations on the linear boundary grid
#include "basis/linear/boundary/operation/datadriven/OperationBLinearBoundary.hpp"
#include "basis/linear/boundary/operation/datadriven/OperationTestLinearBoundary.hpp"
#include "basis/linear/boundary/operation/common/OperationEvalLinearBoundary.hpp"
#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"
// @todo (heinecke) removed this when done
#include "basis/linear/boundary/operation/common/OperationUpDownTestLinearBoundary.hpp"

#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLBLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLDLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLELinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLFLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/OperationLTwoDotProductLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLogLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLogLinearBoundary.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(std::istream& istr) : Grid(istr)
{

}

LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(BoundingBox& BB)
{
	this->storage = new GridStorage(BB);
}

LinearTrapezoidBoundaryGrid::~LinearTrapezoidBoundaryGrid()
{
}

const char* LinearTrapezoidBoundaryGrid::getType()
{
	return "linearTrapezoidBoundary";
}

Grid* LinearTrapezoidBoundaryGrid::unserialize(std::istream& istr)
{
	return new LinearTrapezoidBoundaryGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearTrapezoidBoundaryGrid::createGridGenerator()
{
	return new TrapezoidBoundaryGridGenerator(this->storage);
}

OperationB* LinearTrapezoidBoundaryGrid::createOperationB()
{
	return new OperationBLinearBoundary(this->storage);
}

OperationBVectorized* LinearTrapezoidBoundaryGrid::createOperationBVectorized(const std::string& VecType)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationLaplace()
{
	return new OperationLaplaceLinearBoundary(this->storage);
}

OperationEval* LinearTrapezoidBoundaryGrid::createOperationEval()
{
	return new OperationEvalLinearBoundary(this->storage);
}

OperationTest* LinearTrapezoidBoundaryGrid::createOperationTest()
{
	return new OperationTestLinearBoundary(this->storage);
}

OperationHierarchisation* LinearTrapezoidBoundaryGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearBoundary(this->storage);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationLTwoDotProduct()
{
	return new OperationLTwoDotProductLinearBoundary(this->storage);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationLB()
{
	return new OperationLBLinearBoundary(this->storage);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationLD()
{
	return new OperationLDLinearBoundary(this->storage);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationLE()
{
	return new OperationLELinearBoundary(this->storage);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationLF()
{
	return new OperationLFLinearBoundary(this->storage);
}

// @todo (heinecke) removed this when done
OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationUpDownTest()
{
	return new OperationUpDownTestLinearBoundary(this->storage);
}

// finance operations
/////////////////////
OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationDelta(DataVector& coef)
{
	return new OperationDeltaLinearBoundary(this->storage, coef);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationGamma(DataMatrix& coef)
{
	return new OperationGammaLinearBoundary(this->storage, coef);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationDeltaLog(DataVector& coef)
{
	return new OperationDeltaLogLinearBoundary(this->storage, coef);
}

OperationMatrix* LinearTrapezoidBoundaryGrid::createOperationGammaLog(DataMatrix& coef)
{
	return new OperationGammaLogLinearBoundary(this->storage, coef);
}

}
