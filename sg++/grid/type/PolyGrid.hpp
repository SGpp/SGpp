/*****************************************************************************/
/* This file is part of sg++, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sg++ is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sg++ is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sg++; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

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
	PolyGrid(size_t dim, size_t degree);
	virtual ~PolyGrid();

	virtual const char* getType();
	virtual void serialize(std::ostream& ostr);

	virtual OperationB* createOperationB();
	virtual GridGenerator* createGridGenerator();
	virtual OperationMatrix* createOperationLaplace();
	virtual OperationEval* createOperationEval();
	virtual OperationHierarchisation* createOperationHierarchisation();

	static Grid* unserialize(std::istream& istr);

protected:
	size_t degree;
};

}

#endif /* POLYGRID_HPP */
