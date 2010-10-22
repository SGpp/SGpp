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
#include "basis/linear/noboundary/operation/datadriven/OperationBLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeSSELinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeSPSSELinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeAVXLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeSPAVXLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationTestLinear.hpp"
#include "basis/linear/noboundary/operation/common/OperationEvalLinear.hpp"
#include "basis/linear/noboundary/operation/common/OperationHierarchisationLinear.hpp"

#ifdef USEOCL
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeSPOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeSPHybridSSEOCLLinear.hpp"
#endif

#ifdef USEARBB
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeArBBLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationBIterativeSPArBBLinear.hpp"
#endif

#include "basis/linear/noboundary/operation/pde/OperationLaplaceLinear.hpp"
#include "basis/linear/noboundary/operation/pde/OperationLTwoDotProductLinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLELinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLBLinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLFLinear.hpp"
#include "basis/linear/noboundary/operation/pde/financeHW1D/OperationLDLinear.hpp"

#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLinear.hpp"
#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLinear.hpp"
#include "basis/linear/noboundary/operation/pde/finance/OperationDeltaLogLinear.hpp"
#include "basis/linear/noboundary/operation/pde/finance/OperationGammaLogLinear.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
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

OperationB* LinearGrid::createOperationB()
{
	return new OperationBLinear(this->storage);
}

OperationBVectorized* LinearGrid::createOperationBVectorized(const std::string& VecType)
{
	if (VecType == "SSE")
	{
		return new OperationBIterativeSSELinear(this->storage);
	}
	else if (VecType == "AVX")
	{
		return new OperationBIterativeAVXLinear(this->storage);
	}
#ifdef USEOCL
	else if (VecType == "OCL")
	{
		return new OperationBIterativeOCLLinear(this->storage);
	}
#endif
#ifdef USEARBB
	else if (VecType == "ArBB")
	{
		return new OperationBIterativeArBBLinear(this->storage);
	}
#endif
	else
	{
		throw factory_exception("Unsupported vectorization type");
	}
}

OperationBVectorizedSP* LinearGrid::createOperationBVectorizedSP(const std::string& VecType)
{
	if (VecType == "SSE")
	{
		return new OperationBIterativeSPSSELinear(this->storage);
	}
	else if (VecType == "AVX")
	{
		return new OperationBIterativeSPAVXLinear(this->storage);
	}
#ifdef USEOCL
	else if (VecType == "OCL")
	{
		return new OperationBIterativeSPOCLLinear(this->storage);
	}
	else if (VecType == "HYBRID_SSE_OCL")
	{
		return new OperationBIterativeSPHybridSSEOCLLinear(this->storage);
	}
#endif
#ifdef USEARBB
	else if (VecType == "ArBB")
	{
		return new OperationBIterativeSPArBBLinear(this->storage);
	}
#endif
	else
	{
		throw factory_exception("Unsupported vectorization type");
	}
}

OperationMatrix* LinearGrid::createOperationLaplace()
{
	return new OperationLaplaceLinear(this->storage);
}

OperationEval* LinearGrid::createOperationEval()
{
	return new OperationEvalLinear(this->storage);
}

OperationTest* LinearGrid::createOperationTest()
{
	return new OperationTestLinear(this->storage);
}

OperationHierarchisation* LinearGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinear(this->storage);
}

OperationMatrix* LinearGrid::createOperationLTwoDotProduct()
{
	return new OperationLTwoDotProductLinear(this->storage);
}

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

// @todo (heinecke) removed this when done
OperationMatrix* LinearGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* LinearGrid::createOperationDelta(DataVector& coef)
{
	return new OperationDeltaLinear(this->storage, coef);
}

OperationMatrix* LinearGrid::createOperationGamma(DataMatrix& coef)
{
	return new OperationGammaLinear(this->storage, coef);
}

OperationMatrix* LinearGrid::createOperationDeltaLog(DataVector& coef)
{
	return new OperationDeltaLogLinear(this->storage, coef);
}

OperationMatrix* LinearGrid::createOperationGammaLog(DataMatrix& coef)
{
	return new OperationGammaLogLinear(this->storage, coef);
}

}
