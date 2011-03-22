/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/LinearStretchedGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the linearstretched grid
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSSELinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPSSELinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeAVXLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPAVXLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationTestLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/common/OperationEvalLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/common/OperationHierarchisationLinearStretched.hpp"
//
//#ifdef USEOCL
//#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeOCLLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPOCLLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPHybridSSEOCLLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeHybridSSEOCLLinearStretched.hpp"
//#endif
//
//#ifdef USEARBB
//#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeArBBLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPArBBLinearStretched.hpp"
//#endif

#ifdef USEOCL
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeOCLLinear.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPOCLLinear.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPHybridSSEOCLLinear.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeHybridSSEOCLLinear.hpp"
#endif

#ifdef USEARBB
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeArBBLinear.hpp"
#include "basis/linearstretched/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPArBBLinear.hpp"
#endif

#include "basis/linearstretched/noboundary/operation/pde/OperationLaplaceLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/OperationLTwoDotProductLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLELinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLBLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLFLinearStretched.hpp"
//#include "basis/linearstretched/noboundary/operation/pde/financeHW1D/OperationLDLinearStretched.hpp"

#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationDeltaLogLinearStretched.hpp"
#include "basis/linearstretched/noboundary/operation/pde/finance/OperationGammaLogLinearStretched.hpp"

#include "exception/factory_exception.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

LinearStretchedGrid::LinearStretchedGrid(std::istream& istr) : Grid(istr)
{

}

LinearStretchedGrid::LinearStretchedGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

LinearStretchedGrid::LinearStretchedGrid(Stretching& BB)
{
	this->storage = new GridStorage(BB);
}

LinearStretchedGrid::~LinearStretchedGrid()
{
}

const char* LinearStretchedGrid::getType()
{
	return "linearStretched";
}

Grid* LinearStretchedGrid::unserialize(std::istream& istr)
{
	return new LinearStretchedGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* LinearStretchedGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationMultipleEval* LinearStretchedGrid::createOperationMultipleEval(DataMatrix* dataset)
{
	return new OperationMultipleEvalLinearStretched(this->storage, dataset);
}

OperationMultipleEvalVectorized* LinearStretchedGrid::createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset)
{
	if (VecType == "SSE")
	{
		return new OperationMultipleEvalIterativeSSELinearStretched(this->storage, dataset);
	}
	else if (VecType == "AVX")
	{
		return new OperationMultipleEvalIterativeAVXLinearStretched(this->storage, dataset);
	}
#ifdef USEOCL
	else if (VecType == "OCL")
	{
		return new OperationMultipleEvalIterativeOCLLinearStretched(this->storage, dataset);
	}
	else if (VecType == "HYBRID_SSE_OCL")
	{
		return new OperationMultipleEvalIterativeHybridSSEOCLLinearStretched(this->storage, dataset);
	}
#endif
#ifdef USEARBB
	else if (VecType == "ArBB")
	{
		return new OperationMultipleEvalIterativeArBBLinearStretched(this->storage, dataset);
	}
#endif
	else
	{
		throw factory_exception("Unsupported vectorization type");
	}
}

OperationMultipleEvalVectorizedSP* LinearStretchedGrid::createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset)
{
	if (VecType == "SSE")
	{
		return new OperationMultipleEvalIterativeSPSSELinearStretched(this->storage, dataset);
	}
	else if (VecType == "AVX")
	{
		return new OperationMultipleEvalIterativeSPAVXLinearStretched(this->storage, dataset);
	}
#ifdef USEOCL
	else if (VecType == "OCL")
	{
		return new OperationMultipleEvalIterativeSPOCLLinearStretched(this->storage, dataset);
	}
	else if (VecType == "HYBRID_SSE_OCL")
	{
		return new OperationMultipleEvalIterativeSPHybridSSEOCLLinearStretched(this->storage, dataset);
	}
#endif
#ifdef USEARBB
	else if (VecType == "ArBB")
	{
		return new OperationMultipleEvalIterativeSPArBBLinearStretched(this->storage, dataset);
	}
#endif
	else
	{
		throw factory_exception("Unsupported vectorization type");
	}
}
OperationMatrix* LinearStretchedGrid::createOperationLaplace()
{
	return new OperationLaplaceLinearStretched(this->storage);
}

OperationEval* LinearStretchedGrid::createOperationEval()
{
	return new OperationEvalLinearStretched(this->storage);
}

OperationTest* LinearStretchedGrid::createOperationTest()
{
	return new OperationTestLinearStretched(this->storage);
}

OperationHierarchisation* LinearStretchedGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationLinearStretched(this->storage);
}

OperationMatrix* LinearStretchedGrid::createOperationLTwoDotProduct()
{
	return new OperationLTwoDotProductLinearStretched(this->storage);
}

OperationMatrix* LinearStretchedGrid::createOperationLE()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearStretchedGrid::createOperationLB()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearStretchedGrid::createOperationLF()
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* LinearStretchedGrid::createOperationLD()
{
	throw factory_exception("Unsupported operation");
}

// @todo (heinecke) removed this when done
OperationMatrix* LinearStretchedGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* LinearStretchedGrid::createOperationDelta(DataVector& coef)
{
	return new OperationDeltaLinearStretched(this->storage, coef);
}

OperationMatrix* LinearStretchedGrid::createOperationGamma(DataMatrix& coef)
{
	return new OperationGammaLinearStretched(this->storage, coef);
}

OperationMatrix* LinearStretchedGrid::createOperationDeltaLog(DataVector& coef)
{
	return new OperationDeltaLogLinearStretched(this->storage, coef);
}

OperationMatrix* LinearStretchedGrid::createOperationGammaLog(DataMatrix& coef)
{
	return new OperationGammaLogLinearStretched(this->storage, coef);
}

OperationConvert* LinearStretchedGrid::createOperationConvert()
{
	throw factory_exception("Unsupported operation");
}

}
