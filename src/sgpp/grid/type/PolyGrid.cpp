/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/PolyGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the poly grid
#include "basis/poly/operation/datadriven/OperationMultipleEvalPoly.hpp"
#include "basis/poly/operation/datadriven/OperationTestPoly.hpp"
#include "basis/poly/operation/common/OperationEvalPoly.hpp"
#include "basis/poly/operation/common/OperationHierarchisationPoly.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

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

OperationMultipleEval* PolyGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalPoly(this->storage, this->degree, dataset);
}

OperationMultipleEvalVectorized* PolyGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMultipleEvalVectorizedSP* PolyGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PolyGrid::createOperationLaplace()
{
	throw factory_exception("Unsupported operation");
}

OperationEval* PolyGrid::createOperationEval()
{
	return new OperationEvalPoly(this->storage, this->degree);
}

OperationTest* PolyGrid::createOperationTest()
{
	return new OperationTestPoly(this->storage, this->degree);
}

OperationHierarchisation* PolyGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationPoly(this->storage, this->degree);
}

OperationMatrix* PolyGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PolyGrid::createOperationLB()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PolyGrid::createOperationLD()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PolyGrid::createOperationLE()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PolyGrid::createOperationLF()
{
	throw factory_exception("Unsupported operation");
}

// @todo (heinecke) removed this when done
OperationMatrix* PolyGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* PolyGrid::createOperationDelta(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PolyGrid::createOperationGamma(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PolyGrid::createOperationDeltaLog(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* PolyGrid::createOperationGammaLog(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationConvert* PolyGrid::createOperationConvert()
{
	throw factory_exception("Unsupported operation");
}

}
}
