/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/ModPolyGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the mod poly grid
//#include "basis/modpoly/operation/datadriven/OperationMultipleEvalModPoly.hpp"
//#include "basis/modpoly/operation/datadriven/OperationTestModPoly.hpp"
//#include "basis/modpoly/operation/common/OperationEvalModPoly.hpp"
//#include "basis/modpoly/operation/common/OperationHierarchisationModPoly.hpp"

#include "exception/factory_exception.hpp"


#include <iostream>

namespace sg
{
namespace base
{

ModPolyGrid::ModPolyGrid(std::istream& istr) : Grid(istr), degree(-1)
{
	istr >> degree;
}


ModPolyGrid::ModPolyGrid(size_t dim, size_t degree) : degree(degree)
{
	this->storage = new GridStorage(dim);
}

ModPolyGrid::~ModPolyGrid()
{
}

const char* ModPolyGrid::getType()
{
	return "modpoly";
}

size_t ModPolyGrid::getDegree() const
{
	return this->degree;
}

Grid* ModPolyGrid::unserialize(std::istream& istr)
{
	return new ModPolyGrid(istr);
}

void ModPolyGrid::serialize(std::ostream& ostr)
{
	this->Grid::serialize(ostr);
	ostr << degree << std::endl;
}


/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* ModPolyGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

/*OperationMultipleEval* ModPolyGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalModPoly(this->storage, this->degree, dataset);
}*/

/*OperationMultipleEvalVectorized* ModPolyGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMultipleEvalVectorizedSP* ModPolyGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
{
	throw factory_exception("Unsupported operation");
}*/

/*OperationMatrix* ModPolyGrid::createOperationLaplace()
{
	throw factory_exception("Unsupported operation");
}*/

/*OperationEval* ModPolyGrid::createOperationEval()
{
	return new OperationEvalModPoly(this->storage, this->degree);
}*/

//OperationTest* ModPolyGrid::createOperationTest()
//{
//	return new OperationTestModPoly(this->storage, this->degree);
//}

//OperationHierarchisation* ModPolyGrid::createOperationHierarchisation()
//{
//	return new OperationHierarchisationModPoly(this->storage, this->degree);
//}
/*
OperationMatrix* ModPolyGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}*/
/*
OperationMatrix*  ModPolyGrid::createOperationLB()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModPolyGrid::createOperationLD()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModPolyGrid::createOperationLE()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix*  ModPolyGrid::createOperationLF()
{
	throw factory_exception("Unsupported operation");
}
*/
// @todo (heinecke) removed this when done
/*
OperationMatrix* ModPolyGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}*/

// finance operations
/////////////////////
//OperationMatrix* ModPolyGrid::createOperationDelta(DataVector& coef)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* ModPolyGrid::createOperationGamma(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}*/

//OperationMatrix* ModPolyGrid::createOperationDeltaLog(DataVector& coef)
//{
//	throw factory_exception("Unsupported operation");
//}

/*OperationMatrix* ModPolyGrid::createOperationGammaLog(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}*/

//OperationConvert* ModPolyGrid::createOperationConvert()
//{
//	throw factory_exception("Unsupported operation");
//}

}
}
