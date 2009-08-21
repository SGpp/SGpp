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

#ifndef OPERATIONUPDOWNTESTLINEARBOUNDARYUSCALED_HPP
#define OPERATIONUPDOWNTESTLINEARBOUNDARYUSCALED_HPP

#include "basis/linearboundaryUScaled/algorithm_sweep/PhiPhiDownLinearBoundaryUScaled.hpp"
#include "basis/linearboundaryUScaled/algorithm_sweep/PhiPhiUpLinearBoundaryUScaled.hpp"

#include "operation/OperationMatrix.hpp"

#include "algorithm/UnidirGradient.hpp"
#include "algorithm/sweep.hpp"

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

/**
 * Temporal Test class for Up/Down Algorithm
 */
class OperationUpDownTestLinearBoundaryUScaled: public OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationUpDownTestLinearBoundaryUScaled(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~OperationUpDownTestLinearBoundaryUScaled() {}


	virtual void mult(DataVector& alpha, DataVector& result)
	{
		this->updown(alpha, result);
	}

protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	GridStorage* storage;

	/**
	 * Starting point of the complete up-down scheme
	 *
	 * @param alpha contains the grid points coefficients
	 * @param result contains the result of the laplace operator
	 */
	void updown(DataVector& alpha, DataVector& result)
	{
		DataVector beta(result.getSize());
		result.setAll(0.0);
		beta.setAll(0.0);

		this->updown(alpha, beta, storage->dim() - 1);

		result.add(beta);
	}

	/**
	 * Recursive procedure for updown(). In dimension <i>gradient_dim</i> the L2 scalar product of the
	 * gradients is used. In all other dimensions only the L2 scalar product.
	 *
	 * @param dim the current dimension
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void updown(DataVector& alpha, DataVector& result, size_t dim)
	{
		//Unidirectional scheme
		if(dim > 0)
		{
			// Reordering ups and downs
			// Use previously calculated ups for all future calculations
			// U* -> UU* and UD*

			DataVector temp(alpha.getSize());
			up(alpha, temp, dim);
			updown(temp, result, dim-1);


			// Same from the other direction:
			// *D -> *UD and *DD

			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1);
			down(temp, result_temp, dim);


			//Overall memory use: 2*|alpha|*(d-1)

			result.add(result_temp);
		}
		else
		{
			// Terminates dimension recursion
			up(alpha, result, dim);

			DataVector temp(alpha.getSize());
			down(alpha, temp, dim);

			result.add(temp);
		}
	}

	void up(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::PhiPhiUpLinearBoundaryUScaled func(this->storage);
		sweep<detail::PhiPhiUpLinearBoundaryUScaled> s(func, this->storage);
		s.sweep1D_Boundary(alpha, result, dim);
	}

	void down(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::PhiPhiDownLinearBoundaryUScaled func(this->storage);
		sweep<detail::PhiPhiDownLinearBoundaryUScaled> s(func, this->storage);
		s.sweep1D_Boundary(alpha, result, dim);
	}
};

}

#endif /* OPERATIONUPDOWNTESTLINEARBOUNDARYUSCALED_HPP */
