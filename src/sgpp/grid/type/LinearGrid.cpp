/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/LinearGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the linear grid
//#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalLinear.hpp"

//#include "basis/linear/noboundary/operation/datadriven/OperationTestLinear.hpp"
//#include "basis/linear/noboundary/operation/common/OperationEvalLinear.hpp"
//#include "basis/linear/noboundary/operation/common/OperationHierarchisationLinear.hpp"

//#include "basis/linear/noboundary/operation/pde/OperationLaplaceLinear.hpp"
//#include "basis/linear/noboundary/operation/pde/OperationLTwoDotProductLinear.hpp"
/*#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLELinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLBLinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLFLinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLDLinear.hpp"*/

//#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLinear.hpp"
//#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLinear.hpp"
//#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLogLinear.hpp"
#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLogLinear.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>
using namespace sg::parallel;
using namespace sg::finance;
using namespace sg::pde;
using namespace sg::datadriven;

namespace sg
{
namespace base
{

LinearGrid::LinearGrid(std::istream& istr) : Grid(istr)
{

}

LinearGrid::LinearGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearGrid::LinearGrid(BoundingBox& BB)
{
	this->storage = new GridStorage(BB);
}

LinearGrid::~LinearGrid()
{
}

const char* LinearGrid::getType()
{
	return "linear";
}

Grid* LinearGrid::unserialize(std::istream& istr)
{
	return new LinearGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

/*OperationMultipleEval* LinearGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalLinear(this->storage, dataset);
}*/

/*OperationMultipleEvalVectorized* LinearGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
{
	if (VecType == "SSE")
	{
		return new OperationMultipleEvalIterativeSSELinear(this->storage, dataset);
	}
	else if (VecType == "AVX")
	{
		return new OperationMultipleEvalIterativeAVXLinear(this->storage, dataset);
	}
#ifdef USEOCL
	else if (VecType == "OCL")
	{
		return new OperationMultipleEvalIterativeOCLLinear(this->storage, dataset);
	}
	else if (VecType == "HYBRID_SSE_OCL")
	{
		return new OperationMultipleEvalIterativeHybridSSEOCLLinear(this->storage, dataset);
	}
#endif
#ifdef USEARBB
	else if (VecType == "ArBB")
	{
		return new OperationMultipleEvalIterativeArBBLinear(this->storage, dataset);
	}
#endif
	else
	{
		throw factory_exception("Unsupported vectorization type");
	}
}

OperationMultipleEvalVectorizedSP* LinearGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
{
	if (VecType == "SSE")
	{
		return new OperationMultipleEvalIterativeSPSSELinear(this->storage, dataset);
	}
	else if (VecType == "AVX")
	{
		return new OperationMultipleEvalIterativeSPAVXLinear(this->storage, dataset);
	}
#ifdef USEOCL
	else if (VecType == "OCL")
	{
		return new OperationMultipleEvalIterativeSPOCLLinear(this->storage, dataset);
	}
	else if (VecType == "HYBRID_SSE_OCL")
	{
		return new OperationMultipleEvalIterativeSPHybridSSEOCLLinear(this->storage, dataset);
	}
#endif
#ifdef USEARBB
	else if (VecType == "ArBB")
	{
		return new OperationMultipleEvalIterativeSPArBBLinear(this->storage, dataset);
	}
#endif
	else
	{
		throw factory_exception("Unsupported vectorization type");
	}
}*/

/*
OperationMatrix* LinearGrid::createOperationLaplace()
{
	return new OperationLaplaceLinear(this->storage);
}
*/
/*
OperationEval* LinearGrid::createOperationEval()
{
	return new OperationEvalLinear(this->storage);
}*/

//OperationTest* LinearGrid::createOperationTest()
//{
//	return new OperationTestLinear(this->storage);
//}

//OperationHierarchisation* LinearGrid::createOperationHierarchisation()
//{
//	return new OperationHierarchisationLinear(this->storage);
//}
/*
OperationMatrix* LinearGrid::createOperationLTwoDotProduct()
{
	return new OperationLTwoDotProductLinear(this->storage);
}*/
/*
OperationMatrix* LinearGrid::createOperationLE()
{
	return new OperationLELinear(this->storage);
}

OperationMatrix* LinearGrid::createOperationLB()
{
	return new OperationLBLinear(this->storage);
}

OperationMatrix* LinearGrid::createOperationLF()
{
	return new OperationLFLinear(this->storage);
}

OperationMatrix* LinearGrid::createOperationLD()
{
	return new OperationLDLinear(this->storage);
}
*/

// @todo (heinecke) removed this when done
/*OperationMatrix* LinearGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}*/

// finance operations
/////////////////////
//OperationMatrix* LinearGrid::createOperationDelta(DataVector& coef)
//{
//	return new OperationDeltaLinear(this->storage, coef);
//}

/*OperationMatrix* LinearGrid::createOperationGamma(DataMatrix& coef)
{
	return new OperationGammaLinear(this->storage, coef);
}*/

//OperationMatrix* LinearGrid::createOperationDeltaLog(DataVector& coef)
//{
//	return new OperationDeltaLogLinear(this->storage, coef);
//}

/*OperationMatrix* LinearGrid::createOperationGammaLog(DataMatrix& coef)
{
	return new OperationGammaLogLinear(this->storage, coef);
}*/

//OperationConvert* LinearGrid::createOperationConvert()
//{
//	throw factory_exception("Unsupported operation");
//}

}
}
