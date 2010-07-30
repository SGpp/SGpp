/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef POLYGRID_HPP
#define POLYGRID_HPP

#include "grid/Grid.hpp"

#include <iostream>

namespace sg
{

/**
 * grid with polynomial base functions
 */
class PolyGrid : public Grid
{
protected:
	PolyGrid(std::istream& istr);

public:
	/**
	 * Constructor of grid with polynomial base functions
	 *
	 * @param dim the dimension of the grid
	 * @param degree the max. polynom's degree
	 */
	PolyGrid(size_t dim, size_t degree);

	/**
	 * Destructor
	 */
	virtual ~PolyGrid();

	virtual const char* getType();
	virtual void serialize(std::ostream& ostr);

	virtual OperationB* createOperationB();
	virtual GridGenerator* createGridGenerator();
	virtual OperationMatrix* createOperationLaplace();
	virtual OperationEval* createOperationEval();
	virtual OperationTest* createOperationTest();
	virtual OperationHierarchisation* createOperationHierarchisation();
	virtual OperationMatrix* createOperationLTwoDotProduct();

	// @todo (heinecke) remove this when done
	virtual OperationMatrix* createOperationUpDownTest();

	// finance operations
	virtual OperationMatrix* createOperationDelta(DataVector& coef);
	virtual OperationMatrix* createOperationGamma(DataVector& coef);
	virtual OperationMatrix* createOperationDeltaLog(DataVector& coef);
	virtual OperationMatrix* createOperationGammaLog(DataVector& coef);

	static Grid* unserialize(std::istream& istr);

protected:
	/// max. polynom's degree
	size_t degree;
};

}

#endif /* POLYGRID_HPP */
