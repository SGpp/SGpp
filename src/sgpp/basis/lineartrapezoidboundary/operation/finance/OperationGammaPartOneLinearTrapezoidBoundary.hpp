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

	/**
	 * Recursive procedure for updown().
	 *
	 * @todo (heinecke) add description
	 *
	 * @param dim the current dimension
	 * @param operation_dim_one the dimension in which to use the evaluation
	 * @param operation_dim_two the dimension in which to use the gradient
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown(DataVector& alpha, DataVector& result, size_t dim, size_t operation_dim_one, size_t operation_dim_two);

	/**
	 * All calculations for gradient
	 *
	 * @todo (heinecke) add description
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param operation_dim_one the dimension in which to use the evaluation
	 * @param operation_dim_two the dimension in which to use the gradient
	 */
	void gradient(DataVector& alpha, DataVector& result, size_t dim, size_t operation_dim_one, size_t operation_dim_two);

	/**
	 * All calculations for the function evaluation part of this operation
	 *
	 * @todo (heinecke) add description
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param operation_dim_one the dimension in which to use the evaluation
	 * @param operation_dim_two the dimension in which to use the gradient
	 */
	void SurfaceIntegral(DataVector& alpha, DataVector& result, size_t dim, size_t operation_dim_one, size_t operation_dim_two);

	/**
	 * All calculations for the function evaluation part of this operation, multiplied by squared underlying price
	 *
	 * @todo (heinecke) add description
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param operation_dim_one the dimension in which to use the evaluation
	 * @param operation_dim_two the dimension in which to use the gradient
	 */
	void SurfaceIntegralSquared(DataVector& alpha, DataVector& result, size_t dim, size_t operation_dim_one, size_t operation_dim_two);

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
	 * @todo (heinecke, nice) complete mathematical description
	 *
	 *s @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	void downGradient(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * up-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
	 * in one dimension
	 *
	 * @todo (heinecke, nice) complete mathematical description
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	void upGradient(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * calculate the surface integral in dimension <i>dim</i> multiplied with x
	 *
	 * @todo (heinecke, nice) complete mathematical description
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that the surface integral is calculated
	 */
	void calcSurfaceIntegral(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * calculate the surface integral in dimension <i>dim</i> multiplied with x^2
	 *
	 * @todo (heinecke, nice) complete mathematical description
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that the surface integral is calculated
	 */
	void calcSurfaceIntegralSquared(DataVector& alpha, DataVector& result, size_t dim);
};

}

#endif /* OPERATIONGAMMAPARTONELINEARTRAPEZOIDBOUNDARY_HPP */
