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
//#include "basis/modlinear/operation/datadriven/OperationMultipleEvalModLinear.hpp"
//#include "basis/modlinear/operation/datadriven/OperationTestModLinear.hpp"
//#include "basis/modlinear/operation/common/OperationEvalModLinear.hpp"
//#include "basis/modlinear/operation/common/OperationHierarchisationModLinear.hpp"
//#include "basis/modlinear/operation/pde/OperationLaplaceModLinear.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>
using namespace sg::pde;
using namespace sg::datadriven;

namespace sg
{
namespace base
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

/*OperationMultipleEval* ModLinearGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalModLinear(this->storage, dataset);
}*/

/*OperationMultipleEvalVectorized* ModLinearGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMultipleEvalVectorizedSP* ModLinearGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
{
	throw factory_exception("Unsupported operation");
}*/

/*OperationMatrix* ModLinearGrid::createOperationLaplace()
{
	return new OperationLaplaceModLinear(this->storage);
}*/
/*
OperationEval* ModLinearGrid::createOperationEval()
{
	return new OperationEvalModLinear(this->storage);
}*/

//OperationTest* ModLinearGrid::createOperationTest()
//{
//	return new OperationTestModLinear(this->storage);
//}

//OperationHierarchisation* ModLinearGrid::createOperationHierarchisation()
//{
//	return new OperationHierarchisationModLinear(this->storage);
//}

/*OperationMatrix* ModLinearGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}*/
/*OperationMatrix*  ModLinearGrid::createOperationLB()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModLinearGrid::createOperationLD()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModLinearGrid::createOperationLE()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModLinearGrid::createOperationLF()
{
	throw factory_exception("Unsupported operation");
}*/
// @todo (heinecke) removed this when done
/*
OperationMatrix* ModLinearGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}*/

// finance operations
/////////////////////
//OperationMatrix* ModLinearGrid::createOperationDelta(DataVector& coef)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* ModLinearGrid::createOperationGamma(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}*/

//OperationMatrix* ModLinearGrid::createOperationDeltaLog(DataVector& coef)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* ModLinearGrid::createOperationGammaLog(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}*/

//OperationConvert* ModLinearGrid::createOperationConvert()
//{
//	throw factory_exception("Unsupported operation");
//}

}
}
