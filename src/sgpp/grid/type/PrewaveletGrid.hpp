/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef PREWAVELETGRID_HPP
#define PREWAVELETGRID_HPP

#include "grid/Grid.hpp"

namespace sg
{
namespace base
{

/**
 * grid with prewavelet base functions
 */
class PrewaveletGrid : public Grid
{
protected:
	PrewaveletGrid(std::istream& istr);
	GridStorage* shadowStorage;

public:
	PrewaveletGrid(size_t dim);
	virtual ~PrewaveletGrid();

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

	// finance operations
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

#endif /* PREWAVELETGRID_HPP */
