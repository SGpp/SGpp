/*
 * SquareRootGrid.hpp
 *
 *  Created on: Aug 4, 2010
 *      Author: aliz
 */

#ifndef SQUAREROOTGRID_HPP_
#define SQUAREROOTGRID_HPP_
#include "grid/Grid.hpp"

#include <iostream>
using namespace sg::datadriven;

namespace sg
{
namespace base
{

/**
 * grid with linear base functions with boundaries, pentagon cut
 */
class SquareRootGrid : public Grid
{
protected:
	SquareRootGrid(std::istream& istr);

public:
	/**
	 * Constructor Linear Trapezoid Boundary Grid
	 *
	 * @param dim the dimension of the grid
	 */
	SquareRootGrid(size_t dim);

	/**
	 * Constructor Linear Trapezoid Boundary Grid
	 *
	 * @param BB the BoundingBox of the grid
	 */
	SquareRootGrid(BoundingBox& BB);

	/**
	 * Destructor
	 */
	virtual ~SquareRootGrid();

	virtual const char* getType();

	//virtual OperationMultipleEval* createOperationMultipleEval(DataMatrix* dataset){return 0;};
//	virtual OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset){return 0;};
//	virtual OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(const std::string& VecType, DataMatrixSP* dataset){return 0;};
	virtual GridGenerator* createGridGenerator();
	//virtual OperationMatrix* createOperationLaplace(){return 0;};
//	virtual OperationEval* createOperationEval();
	//virtual OperationTest* createOperationTest(){return 0;};
	//virtual OperationHierarchisation* createOperationHierarchisation();
	//virtual OperationMatrix* createOperationLTwoDotProduct(){return 0;};
	//virtual OperationConvert* createOperationConvert();

	// finance operations
	//virtual OperationMatrix* createOperationDelta(DataVector& coef){return 0;};
	//virtual OperationMatrix* createOperationGamma(DataMatrix& coef){return 0;};
	//virtual OperationMatrix* createOperationDeltaLog(DataVector& coef){return 0;};
	//virtual OperationMatrix* createOperationGammaLog(DataMatrix& coef){return 0;};
	/*virtual OperationMatrix* createOperationLB() {return 0;};
	virtual OperationMatrix* createOperationLD() {return 0;};
	virtual OperationMatrix* createOperationLE() {return 0;};
	virtual OperationMatrix* createOperationLF() {return 0;};*/

	static Grid* unserialize(std::istream& istr);
};

}
}

#endif /* SQUAREROOTGRID_HPP_ */
