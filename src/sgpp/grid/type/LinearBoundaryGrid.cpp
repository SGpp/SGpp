/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"

#include "grid/generation/BoundaryGridGenerator.hpp"

// Include all operations on the linear boundary grid
#include "basis/linear/boundary/operation/datadriven/OperationMultipleEvalLinearBoundary.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSSELinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPSSELinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeAVXLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPAVXLinear.hpp"
#include "basis/linear/boundary/operation/datadriven/OperationTestLinearBoundary.hpp"
#include "basis/linear/boundary/operation/common/OperationEvalLinearBoundary.hpp"
#include "basis/linear/boundary/operation/common/OperationHierarchisationLinearBoundary.hpp"

#ifdef USEOCL
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPHybridSSEOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeHybridSSEOCLLinear.hpp"
#endif

#ifdef USEARBB
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeArBBLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPArBBLinear.hpp"
#endif

#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/OperationLTwoDotProductLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLELinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLBLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLFLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/financeHW1D/OperationLDLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationDeltaLogLinearBoundary.hpp"
#include "basis/linear/boundary/operation/pde/finance/OperationGammaLogLinearBoundary.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{
namespace base
{

LinearBoundaryGrid::LinearBoundaryGrid(std::istream& istr) : Grid(istr)
{

}

LinearBoundaryGrid::LinearBoundaryGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearBoundaryGrid::~LinearBoundaryGrid()
{
}

const char* LinearBoundaryGrid::getType()
{
	return "linearBoundary";
}

Grid* LinearBoundaryGrid::unserialize(std::istream& istr)
{
	return new LinearBoundaryGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearBoundaryGrid::createGridGenerator()
{
	return new BoundaryGridGenerator(this->storage);
}

OperationMultipleEval* LinearBoundaryGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalLinearBoundary(this->storage, dataset);
}

OperationMultipleEvalVectorized* LinearBoundaryGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
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

OperationMultipleEvalVectorizedSP* LinearBoundaryGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
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
}

OperationMatrix* LinearBoundaryGrid::createOperationLaplace()
{
	return new OperationLaplaceLinearBoundary(this->storage);
}

OperationEval* LinearBoundaryGrid::createOperationEval()
{
	return new OperationEvalLinearBoundary(this->storage);
}

OperationTest* LinearBoundaryGrid::createOperationTest()
{
	return new OperationTestLinearBoundary(this->storage);
}

OperationHierarchisation* LinearBoundaryGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearBoundary(this->storage);
}

OperationMatrix* LinearBoundaryGrid::createOperationLTwoDotProduct()
{
	return new OperationLTwoDotProductLinearBoundary(this->storage);
}

OperationMatrix* LinearBoundaryGrid::createOperationLE()
{
	return new OperationLELinearBoundary(this->storage);
}

OperationMatrix* LinearBoundaryGrid::createOperationLB()
{
	return new OperationLBLinearBoundary(this->storage);
}

OperationMatrix* LinearBoundaryGrid::createOperationLF()
{
	return new OperationLFLinearBoundary(this->storage);
}

OperationMatrix* LinearBoundaryGrid::createOperationLD()
{
	return new OperationLDLinearBoundary(this->storage);
}

// @todo (heinecke) removed this when done
OperationMatrix* LinearBoundaryGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* LinearBoundaryGrid::createOperationDelta(DataVector& coef)
{
	return new OperationDeltaLinearBoundary(this->storage, coef);
}

OperationMatrix* LinearBoundaryGrid::createOperationGamma(DataMatrix& coef)
{
	return new OperationGammaLinearBoundary(this->storage, coef);
}

OperationMatrix* LinearBoundaryGrid::createOperationDeltaLog(DataVector& coef)
{
	return new OperationDeltaLogLinearBoundary(this->storage, coef);
}

OperationMatrix* LinearBoundaryGrid::createOperationGammaLog(DataMatrix& coef)
{
	return new OperationGammaLogLinearBoundary(this->storage, coef);
}

OperationConvert* LinearBoundaryGrid::createOperationConvert()
{
	throw factory_exception("Unsupported operation");
}

}
}
