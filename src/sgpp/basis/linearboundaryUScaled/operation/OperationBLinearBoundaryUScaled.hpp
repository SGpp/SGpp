/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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

#ifndef OPERATIONBLINEARBOUNDARYUSCALED_HPP
#define OPERATIONBLINEARBOUNDARYUSCALED_HPP

#include "operation/OperationB.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

/**
 * This class implements OperationB for a grids with linear basis ansatzfunctions
 * with boundaries (pentagon cut through subspace scheme)
 */
class OperationBLinearBoundaryUScaled : public OperationB
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GirdStorage object
	 */
	OperationBLinearBoundaryUScaled(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationBLinearBoundaryUScaled() {}

	virtual void mult(DataVector& alpha, DataVector& data, DataVector& result);
	virtual void multTranspose(DataVector& alpha, DataVector& data, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONBLINEARBOUNDARYUSCALED_HPP */
