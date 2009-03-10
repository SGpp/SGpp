/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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

#ifndef OPERATIONHIERARCHISATIONLINEARBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONLINEARBOUNDARY_HPP

#include "operation/OperationHierarchisation.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * Hierarchisation on sparse grid, linear case
 */
class OperationHierarchisationLinearBoundary : public OperationHierarchisation
{
public:
	OperationHierarchisationLinearBoundary(GridStorage* storage) : storage(storage) {}
	virtual ~OperationHierarchisationLinearBoundary() {}

	virtual void doHierarchisation(DataVector& node_values);
	virtual void doDehierarchisation(DataVector& alpha);

protected:
	GridStorage* storage;
};

}

#endif /* OPERATIONHIERARCHISATIONBOUNDARY_HPP */
