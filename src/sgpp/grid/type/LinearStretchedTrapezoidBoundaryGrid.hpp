/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP
#define LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP

#include "grid/Grid.hpp"

#include <iostream>

namespace sg
{

/**
 * grid with linear base functions with boundaries, pentagon cut
 */
class LinearStretchedTrapezoidBoundaryGrid : public Grid
{
protected:
	LinearStretchedTrapezoidBoundaryGrid(std::istream& istr);

public:
	/**
	 * Constructor Linear Trapezoid Boundary Grid
	 *
	 * @param dim the dimension of the grid
	 */
	LinearStretchedTrapezoidBoundaryGrid(size_t dim);

	/**
	 * Constructor Linear Trapezoid Boundary Grid
	 *
	 * @param BB the Stretching of the grid
	 */
	LinearStretchedTrapezoidBoundaryGrid(Stretching& BB);

	/**
	 * Destructor
	 */
	virtual ~LinearStretchedTrapezoidBoundaryGrid();

	virtual const char* getType();

//	virtual OperationB* createOperationB();
//	virtual OperationBVectorized* createOperationBVectorized(const std::string& VecType);
//	virtual OperationBVectorizedSP* createOperationBVectorizedSP(const std::string& VecType);
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

	// finance operations
	virtual OperationMatrix* createOperationDelta(DataVector& coef);
	virtual OperationMatrix* createOperationGamma(DataMatrix& coef);
	virtual OperationMatrix* createOperationDeltaLog(DataVector& coef);
	virtual OperationMatrix* createOperationGammaLog(DataMatrix& coef);
	// finance operations for hull-white 1D
	virtual OperationMatrix* createOperationLB();
	virtual OperationMatrix* createOperationLD();
	virtual OperationMatrix* createOperationLE();
	virtual OperationMatrix* createOperationLF();


	static Grid* unserialize(std::istream& istr);
};

}

#endif /* LINEARSTRETCHEDTRAPEZOIDBOUNDARYGRID_HPP */
