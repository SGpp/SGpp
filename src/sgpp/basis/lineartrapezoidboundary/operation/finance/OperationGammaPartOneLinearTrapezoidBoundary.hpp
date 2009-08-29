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

#ifndef OPERATIONGAMMAPARTONELINEARTRAPEZOIDBOUNDARY_HPP
#define OPERATIONGAMMAPARTONELINEARTRAPEZOIDBOUNDARY_HPP

#include "operation/common/OperationMatrix.hpp"

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

/**
 * @todo heinecke add description
 */
class OperationGammaPartOneLinearTrapezoidBoundary: public OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param sigma vector that contains the underlyings' standard derivation
	 * @param rho matrix that contains the correlations between the underlyings
	 */
	OperationGammaPartOneLinearTrapezoidBoundary(GridStorage* storage, DataVector& sigma, DataVector rho);

	/**
	 * Destructor
	 */
	virtual ~OperationGammaPartOneLinearTrapezoidBoundary();


	virtual void mult(DataVector& alpha, DataVector& result);

protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	GridStorage* storage;
	/// Pointer to the DataVector of the sigmas
	DataVector* sigmas;
	/// Pointer to the DataVector of the rhos
	DataVector* rhos;

};

}

#endif /* OPERATIONGAMMAPARTONELINEARTRAPEZOIDBOUNDARY_HPP */
