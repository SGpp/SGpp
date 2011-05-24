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
namespace base
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

	//virtual OperationMultipleEval* createOperationMultipleEval(DataMatrix* dataset);
	//virtual OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset);
	//virtual OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset);
	virtual GridGenerator* createGridGenerator();
	//virtual OperationMatrix* createOperationLaplace();
	//virtual OperationEval* createOperationEval();
	//virtual OperationTest* createOperationTest();
	//virtual OperationHierarchisation* createOperationHierarchisation();
	//virtual OperationMatrix* createOperationLTwoDotProduct();
	//virtual OperationConvert* createOperationConvert();

	// finance operations
	//virtual OperationMatrix* createOperationDelta(DataVector& coef);
	//virtual OperationMatrix* createOperationGamma(DataMatrix& coef);
	//virtual OperationMatrix* createOperationDeltaLog(DataVector& coef);
	//virtual OperationMatrix* createOperationGammaLog(DataMatrix& coef);
	// finance operations for hull-white 1D
	/*virtual OperationMatrix* createOperationLB();
	virtual OperationMatrix* createOperationLD();
	virtual OperationMatrix* createOperationLE();
	virtual OperationMatrix* createOperationLF();*/


	static Grid* unserialize(std::istream& istr);
	size_t getDegree() const;

protected:
	/// max. polynom's degree
	size_t degree;
};

}
}

#endif /* POLYGRID_HPP */
