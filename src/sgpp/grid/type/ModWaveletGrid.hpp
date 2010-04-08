/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Dirk Pflueger (pflueged@in.tum.de)                */
/* Copyright (C) 2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

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

};

}

#endif /* MODWAVELETGRID_HPP */
