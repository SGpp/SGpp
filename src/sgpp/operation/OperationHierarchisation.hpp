/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#ifndef OPERATIONHIERARCHISATION_HPP
#define OPERATIONHIERARCHISATION_HPP

#include "data/DataVector.hpp"

namespace sg
{

/**
 * This class implements the hierarchisation and dehierarchisation on the sparse grid
 */
class OperationHierarchisation
{
public:
	/**
	 * Constructor
	 */
	OperationHierarchisation() {}

	/**
	 * Destructor
	 */
	virtual ~OperationHierarchisation() {}

	/**
	 * Implements the hierarchisation on a sprase grid
	 *
	 * @param node_values the functions values in the node base
	 */
	virtual void doHierarchisation(DataVector& node_values) = 0;

	/**
	 * Implements the dehierarchisation on a sprase grid
	 *
	 * @param alpha the coefficients of the sparse grid's base functions
	 */
	virtual void doDehierarchisation(DataVector& alpha) = 0;
};

}

#endif /* OPERATIONHIERARCHISATION_HPP */
