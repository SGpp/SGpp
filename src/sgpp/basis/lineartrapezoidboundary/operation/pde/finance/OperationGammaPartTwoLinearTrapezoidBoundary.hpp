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


#ifndef OPERATIONGAMMAPARTTWOLINEARTRAPEZOIDBOUNDARY_HPP
#define OPERATIONGAMMAPARTTWOLINEARTRAPEZOIDBOUNDARY_HPP

#include "operation/common/OperationMatrix.hpp"

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg
{

/**
 * @todo heinecke add description
 */
class OperationGammaPartTwoLinearTrapezoidBoundary: public OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param sigma vector that contains the underlyings' standard derivation
	 * @param rho matrix that contains the correlations between the underlyings
	 */
	OperationGammaPartTwoLinearTrapezoidBoundary(GridStorage* storage, DataVector& sigma, DataVector& rho);

	/**
	 * Destructor
	 */
	virtual ~OperationGammaPartTwoLinearTrapezoidBoundary();


	virtual void mult(DataVector& alpha, DataVector& result);

private:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	GridStorage* storage;
	/// Pointer to the DataVector of the sigmas
	DataVector* sigmas;
	/// Pointer to the DataVector of the rhos
	DataVector* rhos;

#ifndef USEOMPTHREE
	/**
	 * Recursive procedure for updown().
	 *
	 * @param dim the current dimension
	 * @param gradient_dim the dimension in which to use the gradient
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim);

	/**
	 * All calculations for gradient_dim.
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param gradient_dim the dimension in that the gradient is applied
	 */
	void gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim);
#endif
#ifdef USEOMPTHREE
	/**
	 * Recursive procedure for updown, parallel version using OpenMP 3
	 *
	 * @param dim the current dimension
	 * @param gradient_dim the dimension in which to use the gradient
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim);

	/**
	 * All calculations for gradient_dim, parallel version using OpenMP 3
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param gradient_dim the dimension in that the gradient is applied
	 */
	void gradient_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim);
#endif

	/**
	 * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
	 * Applies the up-part of the one-dimensional mass matrix in one dimension.
	 * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
	 *
	 * @param dim dimension in which to apply the up-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void up(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
	 * Applies the down-part of the one-dimensional mass matrix in one dimension.
	 * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
	 *
	 * @param dim dimension in which to apply the down-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void down(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * down-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	void downGradient(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * up-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	void upGradient(DataVector& alpha, DataVector& result, size_t dim);
};

}

#endif /* OPERATIONGAMMAPARTTWOLINEARTRAPEZOIDBOUNDARY_HPP */
