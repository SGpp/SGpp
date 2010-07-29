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

	virtual OperationB* createOperationB();
	virtual GridGenerator* createGridGenerator();
	virtual OperationEval* createOperationEval();
	virtual OperationTest* createOperationTest();
	virtual OperationHierarchisation* createOperationHierarchisation();

	static Grid* unserialize(std::istream& istr);
};

}

#endif /* LINEARBOUNDARYGRID_HPP */
