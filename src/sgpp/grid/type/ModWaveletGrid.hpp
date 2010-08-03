/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef MODWAVELETGRID_HPP
#define MODWAVELETGRID_HPP

#include "grid/Grid.hpp"

#include <iostream>

namespace sg
{

/**
 * grid with modified polynomial base functions
 */
class ModWaveletGrid : public Grid
{
protected:
	ModWaveletGrid(std::istream& istr);

public:
	/**
	 * Constructor of grid with modified polynomial base functions
	 *
	 * @param dim the dimension of the grid
	 */
	ModWaveletGrid(size_t dim);

	/**
	 * Destructor
	 */
	virtual ~ModWaveletGrid();

	virtual const char* getType();

	virtual OperationB* createOperationB();
	virtual OperationBVectorized* createOperationBVectorized(const std::string& VecType);
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
	virtual OperationMatrix* createOperationGamma(DataMatrix& coef);
	virtual OperationMatrix* createOperationDeltaLog(DataVector& coef);
	virtual OperationMatrix* createOperationGammaLog(DataMatrix& coef);

	static Grid* unserialize(std::istream& istr);

};

}

#endif /* MODWAVELETGRID_HPP */
