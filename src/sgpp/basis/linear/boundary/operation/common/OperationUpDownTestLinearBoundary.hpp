/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
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

#ifndef OPERATIONUPDOWNTESTLINEARTRAPEZOIDBOUNDARY_HPP
#define OPERATIONUPDOWNTESTLINEARTRAPEZOIDBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/SqXdPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/SqXdPhidPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp"

#include "basis/linear/boundary/algorithm_sweep/XPhidPhiDownBBLinearBoundary.hpp"
#include "basis/linear/boundary/algorithm_sweep/XPhidPhiUpBBLinearBoundary.hpp"

#include "operation/common/OperationMatrix.hpp"

#include "algorithm/datadriven/UnidirGradient.hpp"
#include "algorithm/common/sweep.hpp"

namespace sg
{

/**
 * Temporal Test class for Up/Down Algorithm
 *
 * @version $HEAD$
 */
class OperationUpDownTestLinearBoundary: public OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationUpDownTestLinearBoundary(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~OperationUpDownTestLinearBoundary() {}


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
			DataVector temp(alpha.getSize());
			up(alpha, temp, dim);
			updown(temp, result, dim-1);


			// Same from the other direction
			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1);
			down(temp, result_temp, dim);

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
		// phi * phi
		//detail::PhiPhiUpBBLinearBoundary func(this->storage);
		//sweep<detail::PhiPhiUpBBLinearBoundary> s(func, this->storage);

		// x^2 * dphi * dphi
		//detail::SqXdPhidPhiUpBBLinearBoundary func(this->storage);
		//sweep<detail::SqXdPhidPhiUpBBLinearBoundary> s(func, this->storage);

		// x * dphi * phi
		detail::XdPhiPhiUpBBLinearBoundary func(this->storage);
		sweep<detail::XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

		// x * phi * dphi
		//detail::XPhidPhiUpBBLinearBoundary func(this->storage);
		//sweep<detail::XPhidPhiUpBBLinearBoundary> s(func, this->storage);

		s.sweep1D_Boundary(alpha, result, dim);
	}

	void down(DataVector& alpha, DataVector& result, size_t dim)
	{
		// phi * phi
		//detail::PhiPhiDownBBLinearBoundary func(this->storage);
		//sweep<detail::PhiPhiDownBBLinearBoundary> s(func, this->storage);

		// x^2 * dphi * dphi
		//detail::SqXdPhidPhiDownBBLinearBoundary func(this->storage);
		//sweep<detail::SqXdPhidPhiDownBBLinearBoundary> s(func, this->storage);

		// x * dphi * phi
		detail::XdPhiPhiDownBBLinearBoundary func(this->storage);
		sweep<detail::XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

		// x * phi * dphi
		//detail::XPhidPhiDownBBLinearBoundary func(this->storage);
		//sweep<detail::XPhidPhiDownBBLinearBoundary> s(func, this->storage);

		s.sweep1D_Boundary(alpha, result, dim);
	}
};

}

#endif /* OPERATIONUPDOWNTESTLINEARBOUNDARY_HPP */
