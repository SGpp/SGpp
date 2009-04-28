/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007 Jörg Blank (blankj@in.tum.de)                          */
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

#ifndef OPERATIONLAPLACELINEAR_HPP
#define OPERATIONLAPLACELINEAR_HPP

#include "basis/linear/algorithm_sweep/LaplaceDownLinear.hpp"
#include "basis/linear/algorithm_sweep/LaplaceUpLinear.hpp"

#include "operation/OperationMatrix.hpp"

#include "algorithm/UnidirGradient.hpp"
#include "algorithm/sweep.hpp"

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

namespace sg
{

/**
 * Implementation for linear functions of Laplace Operation, linear grids without boundaries
 */
class OperationLaplaceLinear: public OperationMatrix, public UnidirGradient
{
public:
	/**
	 * Construtor of OperationLaplaceLinear
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationLaplaceLinear(GridStorage* storage) : UnidirGradient(storage)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~OperationLaplaceLinear() {}

	virtual void mult(DataVector& alpha, DataVector& result)
	{
		this->updown(alpha, result);
	}

protected:
	virtual void gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
	{
		// In direction gradient_dim we only calculate the norm of the gradient
		// The up-part is empty, thus omitted
		if(dim > 0)
		{
			DataVector temp(alpha.getSize());
			updown(alpha, temp, dim-1, gradient_dim);
			downGradient(temp, result, gradient_dim);
		}
		else
		{
			// Terminates dimension recursion
			downGradient(alpha, result, gradient_dim);
		}
	}

	virtual void up(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::LaplaceUpLinear func(this->storage);
		sweep<detail::LaplaceUpLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}

	virtual void down(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::LaplaceDownLinear func(this->storage);
		sweep<detail::LaplaceDownLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}

	virtual void downGradient(DataVector& alpha, DataVector& result, size_t dim)
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			//only affects the diagonal of the stiffness matrix
			result[i] = alpha[i]*pow(2.0, level+1);
		}
	}

	virtual void upGradient(DataVector& alpha, DataVector& result, size_t dim) {}

};

}

#endif /* OPERATIONLAPLACELINEAR_HPP */
