/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LINEARTRAPEZOIDBOUNDARYGRID_HPP
#define LINEARTRAPEZOIDBOUNDARYGRID_HPP

#include "grid/Grid.hpp"

#include <iostream>

namespace sg
{

/**
 * grid with linear base functions with boundaries, pentagon cut
 */
class LinearTrapezoidBoundaryGrid : public Grid
{
protected:
	LinearTrapezoidBoundaryGrid(std::istream& istr);

public:
	/**
	 * Constructor Linear Trapezoid Boundary Grid
	 *
	 * @param dim the dimension of the grid
	 */
	LinearTrapezoidBoundaryGrid(size_t dim);

	/**
	 * Constructor Linear Trapezoid Boundary Grid
	 *
	 * @param BB the BoundingBox of the grid
	 */
	LinearTrapezoidBoundaryGrid(BoundingBox& BB);

	/**
	 * Destructor
	 */
	virtual ~LinearTrapezoidBoundaryGrid();

	virtual const char* getType();

	virtual OperationB* createOperationB();
	virtual OperationBVectorized* createOperationBVectorized(const std::string& VecType);
	virtual GridGenerator* createGridGenerator();
	virtual OperationMatrix* createOperationLaplace();
	virtual OperationEval* createOperationEval();
	virtual OperationTest* createOperationTest();
	virtual OperationHierarchisation* createOperationHierarchisation();
	virtual OperationMatrix* createOperationLTwoDotProduct();
	virtual OperationMatrix* createOperationLB();
	virtual OperationMatrix* createOperationLD();
	virtual OperationMatrix* createOperationLE();
	virtual OperationMatrix* createOperationLF();

	// @todo (heinecke) remove this when done
	virtual OperationMatrix* createOperationUpDownTest();

	// finance operations
	virtual OperationMatrix* createOperationDelta(DataVector& coef);
	virtual OperationMatrix* createOperationGamma(DataMatrix& coef);
	virtual OperationMatrix* createOperationDeltaLog(DataVector& coef);
	virtual OperationMatrix* createOperationGammaLog(DataMatrix& coef);

	static Grid* unserialize(std::istream& istr);
};

}

#endif /* LINEARTRAPEZOIDBOUNDARYGRID_HPP */
