/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2010 Dirk Pflueger (pflueged@in.tum.de)                     */
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

#ifndef MODBSPLINEGRID_HPP
#define MODBSPLINEGRID_HPP

#include "grid/Grid.hpp"

#include <iostream>

namespace sg
{

/**
 * Grid with modified Bspline basis functions
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
