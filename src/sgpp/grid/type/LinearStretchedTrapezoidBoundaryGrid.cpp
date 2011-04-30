/******************************************************************************
 * Copyright (C) 2010 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author

#include "grid/Grid.hpp"
#include "grid/type/LinearStretchedTrapezoidBoundaryGrid.hpp"

#include "grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
//
//// Include all operations on the linear boundary grid
#include "basis/linearstretched/boundary/operation/datadriven/OperationMultipleEvalLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/datadriven/OperationTestLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/common/OperationEvalLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/common/OperationHierarchisationLinearStretchedBoundary.hpp"
// @todo (heinecke) removed this when done
#include "basis/linearstretched/boundary/operation/common/OperationUpDownTestLinearStretchedBoundary.hpp"

#include "basis/linearstretched/boundary/operation/pde/OperationLaplaceLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLBLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLDLinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLELinearStretchedBoundary.hpp"
//#include "basis/linearstretched/boundary/operation/pde/financeHW1D/OperationLFLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/OperationLTwoDotProductLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationDeltaLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationGammaLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationDeltaLogLinearStretchedBoundary.hpp"
#include "basis/linearstretched/boundary/operation/pde/finance/OperationGammaLogLinearStretchedBoundary.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{
namespace base
{

LinearStretchedTrapezoidBoundaryGrid::LinearStretchedTrapezoidBoundaryGrid(std::istream& istr) : Grid(istr)
{

}

LinearStretchedTrapezoidBoundaryGrid::LinearStretchedTrapezoidBoundaryGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearStretchedTrapezoidBoundaryGrid::LinearStretchedTrapezoidBoundaryGrid(Stretching& BB)
{
//	std::cout<<"creating new grid storage\n";
	this->storage = new GridStorage(BB);
}

LinearStretchedTrapezoidBoundaryGrid::~LinearStretchedTrapezoidBoundaryGrid()
{
}

const char* LinearStretchedTrapezoidBoundaryGrid::getType()
{
	return "linearStretchedTrapezoidBoundary";
}

Grid* LinearStretchedTrapezoidBoundaryGrid::unserialize(std::istream& istr)
{
	return new LinearStretchedTrapezoidBoundaryGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearStretchedTrapezoidBoundaryGrid::createGridGenerator()
{
	return new StretchedTrapezoidBoundaryGridGenerator(this->storage);
}

OperationMultipleEval* LinearStretchedTrapezoidBoundaryGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalLinearStretchedBoundary(this->storage, dataset);
}

OperationMultipleEvalVectorized* LinearStretchedTrapezoidBoundaryGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMultipleEvalVectorizedSP* LinearStretchedTrapezoidBoundaryGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationLaplace()
{
	return new OperationLaplaceLinearStretchedBoundary(this->storage);
}

OperationEval* LinearStretchedTrapezoidBoundaryGrid::createOperationEval()
{
	return new OperationEvalLinearStretchedBoundary(this->storage);
}

OperationTest* LinearStretchedTrapezoidBoundaryGrid::createOperationTest()
{
	return new OperationTestLinearStretchedBoundary(this->storage);
}

OperationHierarchisation* LinearStretchedTrapezoidBoundaryGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearStretchedBoundary(this->storage);
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationLTwoDotProduct()
{
	return new OperationLTwoDotProductLinearStretchedBoundary(this->storage);
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationLB()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationLD()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationLE()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationLF()
{
	throw factory_exception("Unsupported operation");
}

// @todo (heinecke) removed this when done
OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationUpDownTest()
{
	return new OperationUpDownTestLinearStretchedBoundary(this->storage);
}

// finance operations
/////////////////////
OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationDelta(DataVector& coef)
{
	return new OperationDeltaLinearStretchedBoundary(this->storage, coef);
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationGamma(DataMatrix& coef)
{
	return new OperationGammaLinearStretchedBoundary(this->storage, coef);
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationDeltaLog(DataVector& coef)
{
	return new OperationDeltaLogLinearStretchedBoundary(this->storage, coef);
}

OperationMatrix* LinearStretchedTrapezoidBoundaryGrid::createOperationGammaLog(DataMatrix& coef)
{
	return new OperationGammaLogLinearStretchedBoundary(this->storage, coef);
}

OperationConvert* LinearStretchedTrapezoidBoundaryGrid::createOperationConvert()
{
	throw factory_exception("Unsupported operation");
}

}
}
