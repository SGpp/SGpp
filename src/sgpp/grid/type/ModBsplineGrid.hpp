/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef MODBSPLINEGRID_HPP
#define MODBSPLINEGRID_HPP

#include "grid/Grid.hpp"

#include <iostream>

namespace sg
{

/**
 * Grid with modified Bspline basis functions
 * @todo (pflueged) include for factory exception is missing in several classes which use it. It only works, as it is include by a header loaded previously.
 */
class ModBsplineGrid : public Grid
{
protected:
	ModBsplineGrid(std::istream& istr);

public:
	/**
	 * Constructor of grid with modified bspline basis functions
	 *
	 * @param dim the dimension of the grid
     * @param degree the bspline's degree
	 */
	ModBsplineGrid(size_t dim, size_t degree);

	/**
	 * Destructor
	 */
	virtual ~ModBsplineGrid();

	virtual const char* getType();

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

	static Grid* unserialize(std::istream& istr);

	virtual void serialize(std::ostream& ostr);

protected:
    // degree of Bspline
    size_t degree;

};

}

#endif /* MODBSPLINEGRID_HPP */
