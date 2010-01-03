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

#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/PhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/SqXdPhidPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/SqXdPhidPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/XPhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/XPhiPhiUpBBLinearTrapezoidBoundary.hpp"

#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiDownBBLinearTrapezoidBoundary.hpp"
#include "basis/lineartrapezoidboundary/algorithm_sweep/XdPhiPhiUpBBLinearTrapezoidBoundary.hpp"

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
class OperationUpDownTestLinearTrapezoidBoundary: public OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationUpDownTestLinearTrapezoidBoundary(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~OperationUpDownTestLinearTrapezoidBoundary() {}


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
		//detail::PhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
		//sweep<detail::PhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

		// x^2 * dphi * dphi
		//detail::SqXdPhidPhiUpBBLinearTrapezoidBoundary func(this->storage);
		//sweep<detail::SqXdPhidPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

		// x * phi * phi
		//detail::XPhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
		//sweep<detail::XPhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

		// x * dphi * phi
		detail::XdPhiPhiUpBBLinearTrapezoidBoundary func(this->storage);
		sweep<detail::XdPhiPhiUpBBLinearTrapezoidBoundary> s(func, this->storage);

		s.sweep1D_Boundary(alpha, result, dim);
	}

	void down(DataVector& alpha, DataVector& result, size_t dim)
	{
		// phi * phi
		//detail::PhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
		//sweep<detail::PhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

		// x^2 * dphi * dphi
		//detail::SqXdPhidPhiDownBBLinearTrapezoidBoundary func(this->storage);
		//sweep<detail::SqXdPhidPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

		// x * phi * phi
		//detail::XPhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
		//sweep<detail::XPhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

		// x * dphi * phi
		detail::XdPhiPhiDownBBLinearTrapezoidBoundary func(this->storage);
		sweep<detail::XdPhiPhiDownBBLinearTrapezoidBoundary> s(func, this->storage);

		s.sweep1D_Boundary(alpha, result, dim);
	}
};

}

#endif /* OPERATIONUPDOWNTESTLINEARTRAPEZOIDBOUNDARY_HPP */
