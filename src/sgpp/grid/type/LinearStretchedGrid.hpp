/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef LINEARSTRETCHEDGRID_HPP
#define LINEARSTRETCHEDGRID_HPP

#include "grid/Grid.hpp"
//#include "grid/common/BoundingBox.hpp"
#include "grid/common/Stretching.hpp"

#include <iostream>

namespace sg
{
namespace base
{

/**
 * grid with linearstretched base functions
 */
class LinearStretchedGrid : public Grid
{
protected:
	LinearStretchedGrid(std::istream& istr);

public:
	/**
	 * Constructor LinearStretched Grid without boundaries
	 *
	 * @param dim the dimension of the grid
	 */
	LinearStretchedGrid(size_t dim);

	/**
	 * Constructor LinearStretched Grid
	 *
	 * @param BB the BoundingBox of the grid
	 */
	LinearStretchedGrid(Stretching& BB);

	/**
	 * Destructor
	 */
	virtual ~LinearStretchedGrid();

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
	virtual OperationMatrix* createOperationLE();
	virtual OperationMatrix* createOperationLB();
	virtual OperationMatrix* createOperationLF();
	virtual OperationMatrix* createOperationLD();

	static Grid* unserialize(std::istream& istr);
};

}
}

#endif /* LINEARSTRETCHEDGRID_HPP */
