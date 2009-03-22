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

#ifndef OPERATIONHIERARCHISATIONPOLY_HPP
#define OPERATIONHIERARCHISATIONPOLY_HPP

#include "operation/OperationHierarchisation.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

/**
 * Hierarchisation on sparse grid, poly case
 */
class OperationHierarchisationPoly : public OperationHierarchisation
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the polynom's max. degree
	 */
	OperationHierarchisationPoly(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationHierarchisationPoly() {}

	/**
	 * Implements the hierarchisation on a sprase grid with poly base functions
	 *
	 * @param node_values the functions values in the node base
	 *
	 * @todo Implement the hierarchisation on the sparse grid with poly base functions
	 */
	virtual void doHierarchisation(DataVector& node_values);

	/**
	 * Implements the dehierarchisation on a sprase grid with poly base functions
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 *
	 * @todo Implement the dehierarchisation on the sparse grid with poly base functions
	 */
	virtual void doDehierarchisation(DataVector& alpha);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Poly Basis object
	SPolyBase base;
};

}

#endif /* OPERATIONHIERARCHISATIONPOLY_HPP */
