/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef UPDOWNTWOOPDIMS_HPP
#define UPDOWNTWOOPDIMS_HPP

#include "grid/GridStorage.hpp"

#include "operation/common/OperationMatrix.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg
{

/**
 * Implements the Up/Down scheme with two dimensions with special operations: i,j
 *
 * Parallelization with OpenMP 2 / 3 is supported!
 *
 * Only symmetric operations are support --> only
 * the "left lower triangular matrix", i <= j, is calculated, please
 * keep that in mind when designing the coefficient vector:
 * the non-diagonal elements must be multiplied by 2
 * before executing this Up/down scheme!
 *
 * @version $HEAD$
 */
class UpDownTwoOpDims: public OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param coef vector that contains the constant coefficients of this operation
	 */
	UpDownTwoOpDims(GridStorage* storage, DataMatrix& coef);

	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	UpDownTwoOpDims(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~UpDownTwoOpDims();


	virtual void mult(DataVector& alpha, DataVector& result);

protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	GridStorage* storage;
	/// Pointer to the coefficients of this bilinear form
	DataMatrix* coefs;

#ifndef USEOMPTHREE
	/**
	 * Recursive procedure for updown
	 *
	 * @param dim the current dimension
	 * @param op_dim_one the dimension in which to use the first gradient
	 * @param op_dim_two the dimension in which to use the second gradient
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);

	/**
	 * All calculations for gradient
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param op_dim_one the dimension in which to use the first gradient
	 * @param op_dim_two the dimension in which to use the second gradient
	 */
	void specialOpOne(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);

	/**
	 * All calculations for gradient, Part 2
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param op_dim_one the dimension in which to use the first gradient
	 * @param op_dim_two the dimension in which to use the second gradient
	 */
	void specialOpTwo(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);

	/**
	 * if the current dimension is equal to the both special operation dimensions
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param op_dim_one the dimension in which to use the first gradient
	 * @param op_dim_two the dimension in which to use the second gradient
	 */
	void specialOpOneAndOpTwo(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);
#endif
#ifdef USEOMPTHREE
	/**
	 * Recursive procedure for updown, parallel version using OpenMP 3
	 *
	 * @param dim the current dimension
	 * @param op_dim_one the dimension in which to use the first gradient
	 * @param op_dim_two the dimension in which to use the second gradient
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);

	/**
	 * All calculations for gradient, parallel version using OpenMP 3
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param op_dim_one the dimension in which to use the first gradient
	 * @param op_dim_two the dimension in which to use the second gradient
	 */
	void specialOpOne_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);

	/**
	 * All calculations for gradient, Part 2, parallel version using OpenMP 3
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param op_dim_one the dimension in which to use the first gradient
	 * @param op_dim_two the dimension in which to use the second gradient
	 */
	void specialOpTwo_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);

	/**
	 * if the current dimension is equal to the both special operation dimensions
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param op_dim_one the dimension in which to use the first gradient
	 * @param op_dim_two the dimension in which to use the second gradient
	 */
	void specialOpOneAndOpTwo_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);
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
	virtual void up(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
	 * Applies the down-part of the one-dimensional mass matrix in one dimension.
	 * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
	 *
	 * @param dim dimension in which to apply the down-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void down(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * 1D down if the current dim is equal to i
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downOpDimOne(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * 1D up if the current dim is equal to i
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upOpDimOne(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * 1D down if the current dim is equal to j
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * 1D up if the current dim is equal to j
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * 1D down, if the current dim is equal to i and j
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * 1D up, if the current dim is equal to i and j
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim) = 0;
};

}

#endif /* UPDOWNTWOOPDIMS_HPP */
