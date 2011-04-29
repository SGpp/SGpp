/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "grid/Grid.hpp"
#include "grid/type/ModBsplineGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the mod bspline grid
//#include "basis/modbspline/operation/datadriven/OperationMultipleEvalModBspline.hpp"
//#include "basis/modbspline/operation/datadriven/OperationTestModBspline.hpp"
//#include "basis/modbspline/operation/common/OperationEvalModBspline.hpp"
//#include "basis/modbspline/operation/common/OperationHierarchisationModBspline.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{
namespace base
{

ModBsplineGrid::ModBsplineGrid(std::istream& istr) : Grid(istr), degree(-1)
{
    istr >> degree;
}


ModBsplineGrid::ModBsplineGrid(size_t dim, size_t degree) : degree(degree)
{
    this->storage = new GridStorage(dim);
}

ModBsplineGrid::~ModBsplineGrid()
{
}

const char* ModBsplineGrid::getType()
{
	return "modBspline";
}

size_t ModBsplineGrid::getDegree()
{
	return this->degree;
}

Grid* ModBsplineGrid::unserialize(std::istream& istr)
{
	return new ModBsplineGrid(istr);
}

void ModBsplineGrid::serialize(std::ostream& ostr)
{
    this->Grid::serialize(ostr);
    ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* ModBsplineGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

/*OperationMultipleEval* ModBsplineGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalModBspline(this->storage, this->degree, dataset);
}*/

/*OperationMultipleEvalVectorized* ModBsplineGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMultipleEvalVectorizedSP* ModBsplineGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
{
	throw factory_exception("Unsupported operation");
}'/

/*OperationMatrix* ModBsplineGrid::createOperationLaplace()
{
	throw factory_exception("Unsupported operation");
}*/

/*OperationEval* ModBsplineGrid::createOperationEval()
{
	return new OperationEvalModBspline(this->storage, this->degree);
}*/

//OperationTest* ModBsplineGrid::createOperationTest()
//{
//    return new OperationTestModBspline(this->storage, this->degree);
//}

//OperationHierarchisation* ModBsplineGrid::createOperationHierarchisation()
//{
//	return new OperationHierarchisationModBspline(this->storage, this->degree);
//}
/*
OperationMatrix* ModBsplineGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}*/
/*OperationMatrix*  ModBsplineGrid::createOperationLB()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModBsplineGrid::createOperationLD()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModBsplineGrid::createOperationLE()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModBsplineGrid::createOperationLF()
{
	throw factory_exception("Unsupported operation");
}*/
/*
OperationMatrix* ModBsplineGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}*/

// finance operations
/////////////////////
//OperationMatrix* ModBsplineGrid::createOperationDelta(DataVector& coef)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* ModBsplineGrid::createOperationGamma(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}*/

//OperationMatrix* ModBsplineGrid::createOperationDeltaLog(DataVector& coef)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* ModBsplineGrid::createOperationGammaLog(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}*/
//
//OperationConvert* ModBsplineGrid::createOperationConvert()
//{
//	throw factory_exception("Unsupported operation");
//}

}
}
