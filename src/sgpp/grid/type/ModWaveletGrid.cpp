/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/ModWaveletGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the mod wavelet grid
#include "basis/modwavelet/operation/datadriven/OperationMultipleEvalModWavelet.hpp"
#include "basis/modwavelet/operation/datadriven/OperationTestModWavelet.hpp"
#include "basis/modwavelet/operation/common/OperationEvalModWavelet.hpp"
#include "basis/modwavelet/operation/common/OperationHierarchisationModWavelet.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{
namespace base
{

ModWaveletGrid::ModWaveletGrid(std::istream& istr) : Grid(istr)
{
}

ModWaveletGrid::ModWaveletGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

ModWaveletGrid::~ModWaveletGrid()
{
}

const char* ModWaveletGrid::getType()
{
	return "modWavelet";
}

Grid* ModWaveletGrid::unserialize(std::istream& istr)
{
	return new ModWaveletGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* ModWaveletGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationMultipleEval* ModWaveletGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalModWavelet(this->storage, dataset);
}

OperationMultipleEvalVectorized* ModWaveletGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMultipleEvalVectorizedSP* ModWaveletGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationLaplace()
{
	throw factory_exception("Unsupported operation");
}

OperationEval* ModWaveletGrid::createOperationEval()
{
	return new OperationEvalModWavelet(this->storage);
}

OperationTest* ModWaveletGrid::createOperationTest()
{
	return new OperationTestModWavelet(this->storage);
}

OperationHierarchisation* ModWaveletGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationModWavelet(this->storage);
}

OperationMatrix* ModWaveletGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationLB()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationLD()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationLE()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationLF()
{
	throw factory_exception("Unsupported operation");
}

// @todo (heinecke) removed this when done
OperationMatrix* ModWaveletGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* ModWaveletGrid::createOperationDelta(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationGamma(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationDeltaLog(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationGammaLog(DataMatrix& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationConvert* ModWaveletGrid::createOperationConvert()
{
	throw factory_exception("Unsupported operation");
}


}
}
