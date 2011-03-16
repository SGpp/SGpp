/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LINEARBOUNDARYGRID_HPP
#define LINEARBOUNDARYGRID_HPP

#include "grid/Grid.hpp"

#include <iostream>

namespace sg
{
namespace base
{

/**
 * grid with linear base functions with boundaries
 */
class LinearBoundaryGrid : public Grid
{
protected:
	LinearBoundaryGrid(std::istream& istr);

public:
	/**
	 * Constructor for the Linear Boundary Grid
	 *
	 * @param dim the dimension of the grid
	 */
	LinearBoundaryGrid(size_t dim);

	/**
	 * Destructor
	 */
	virtual ~LinearBoundaryGrid();

	virtual const char* getType();

	virtual OperationMultipleEval* createOperationMultipleEval(DataMatrix* dataset);
	virtual OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset);
	virtual OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset);
	virtual GridGenerator* createGridGenerator();
	virtual OperationMatrix* createOperationLaplace();
	virtual OperationEval* createOperationEval();
	virtual OperationTest* createOperationTest();
	virtual OperationHierarchisation* createOperationHierarchisation();
	virtual OperationMatrix* createOperationLTwoDotProduct();
	virtual OperationConvert* createOperationConvert();


	// @todo (heinecke) remove this when done
	virtual OperationMatrix* createOperationUpDownTest();

	// finance operations for Black-Scholes nD
	virtual OperationMatrix* createOperationDelta(DataVector& coef);
	virtual OperationMatrix* createOperationGamma(DataMatrix& coef);
	virtual OperationMatrix* createOperationDeltaLog(DataVector& coef);
	virtual OperationMatrix* createOperationGammaLog(DataMatrix& coef);
	// finance operations for Hull-White 1D
	virtual OperationMatrix* createOperationLE();
	virtual OperationMatrix* createOperationLB();
	virtual OperationMatrix* createOperationLF();
	virtual OperationMatrix* createOperationLD();

	static Grid* unserialize(std::istream& istr);
};

}
}

#endif /* LINEARBOUNDARYGRID_HPP */
