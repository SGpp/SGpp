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

#ifndef OPERATIONHIERARCHISATIONMODBSPLINE_HPP
#define OPERATIONHIERARCHISATIONMODBSPLINE_HPP

#include "operation/common/OperationHierarchisation.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

/**
 * Hierarchisation on sparse grid, mod bspline case
 *
 * @version $HEAD$
 */
class OperationHierarchisationModBspline : public OperationHierarchisation
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the bsplinenom's max. degree
	 */
	OperationHierarchisationModBspline(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationHierarchisationModBspline() {}

	/**
	 * Implements the hierarchisation on a sprase grid with mod bspline base functions
	 *
	 * @param node_values the functions values in the node base
	 *
	 * @todo Implement the hierarchisation on the sparse grid with mod bspline base functions
	 */
	virtual void doHierarchisation(DataVector& node_values);

	/**
	 * Implements the dehierarchisation on a sprase grid with mod bspline base functions
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 *
	 * @todo Implement the dehierarchisation on the sparse grid with mod bspline base functions
	 */
	virtual void doDehierarchisation(DataVector& alpha);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Mod Bspline Basis object
	SModBsplineBase base;
};

}

#endif /* OPERATIONHIERARCHISATIONMODBSPLINE_HPP */
