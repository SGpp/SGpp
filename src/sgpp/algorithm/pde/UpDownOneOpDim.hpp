/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef UPDOWNONEOPDIM_HPP
#define UPDOWNONEOPDIM_HPP

#include <vector>

#include "grid/GridStorage.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"

namespace sg
{

/**
 * Implements the Up/Down scheme with one dimension with a special operation
 *
 * Parallelization with OpenMP 2 / 3 is supported!
 *
 * @version $HEAD$
 */
class UpDownOneOpDim: public OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param coef reference to a DataVector object that contains the bilinear form's constant coefficients; one per dimension
	 */
	UpDownOneOpDim(GridStorage* storage, DataVector& coef);

	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	UpDownOneOpDim(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~UpDownOneOpDim();

	virtual void mult(DataVector& alpha, DataVector& result);

	/**
	 * this functions provides the same functionality as the normal mult routine.
	 * However, it doesn't set up an OpenMP task initialization as the mult routine.
	 * This method has to be called within a OpenMP task parallelized region.
	 *
	 * Using this function is useful in following case: Assuming the solver of a certain
	 * requires several operators in the space discretization (e.g. Black Scholes Equations)
	 * this method can be used to parallelize their calculation which might results results
	 * in a better parallel efficiency on systems with 4 or more cores hence fewer barriers
	 * are needed.
	 *
	 * For full calculation in of mult serval number of up/downs are needed. This number
	 * is equal to the number of the grid's dimensions. All different steps can be executed
	 * in parallel, so here only one up/Down is executed, identified by its special dimension.
	 *
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 * @param operationDim Dimension in which the special operator is applied
	 */
	void multParallelBuildingBlock(DataVector& alpha, DataVector& result, size_t operationDim);


protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	GridStorage* storage;
	/// Pointer to the DataVector of the coefs
	DataVector* coefs;
	/// algorithmic dimensions, operator is applied in this dimensions
	const std::vector<size_t> algoDims;
	/// number of algorithmic dimensions
	const size_t numAlgoDims_;
	/// max number of parallel stages (dimension recursive calls)
	static const size_t maxParallelDims_ = 2;

	/**
	 * Recursive procedure for updown(), parallel version using OpenMP 3
	 *
	 * @param dim the current dimension
	 * @param op_dim the dimension in which a special operation is applied
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim);

	/**
	 * All calculations for gradient_dim, parallel version using OpenMP 3
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param op_dim the dimension in that a special operation is applied
	 */
	virtual void specialOP(DataVector& alpha, DataVector& result, size_t dim, size_t op_dim);

	/**
	 * std 1D up operation
	 *
	 * @param dim dimension in which to apply the up-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void up(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * std 1D down operation
	 *
	 * @param dim dimension in which to apply the down-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void down(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * special 1D down operation that is only executed in one direction
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downOpDim(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * special 1D up operation that is only executed in one direction
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upOpDim(DataVector& alpha, DataVector& result, size_t dim) = 0;
};

}

#endif /* UPDOWNONEOPDIM_HPP */
