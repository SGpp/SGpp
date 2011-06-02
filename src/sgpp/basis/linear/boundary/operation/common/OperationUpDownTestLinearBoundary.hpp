/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

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

#include "algorithm/common/sweep.hpp"

namespace sg
{
namespace pde
{

/**
 * Test class for Up/Down Algorithms
 *
 * @version $HEAD$
 */
class OperationUpDownTestLinearBoundary: public sg::base::OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's sg::base::GridStorage object
	 */
	OperationUpDownTestLinearBoundary(sg::base::GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~OperationUpDownTestLinearBoundary();

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

protected:
	typedef sg::base::GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	sg::base::GridStorage* storage;

	/**
	 * Starting point of the complete up-down scheme
	 *
	 * @param alpha contains the grid points coefficients
	 * @param result contains the result of the laplace operator
	 */
	void updown(sg::base::DataVector& alpha, sg::base::DataVector& result);

	/**
	 * Recursive procedure for updown(). In dimension <i>gradient_dim</i> the L2 scalar product of the
	 * gradients is used. In all other dimensions only the L2 scalar product.
	 *
	 * @param dim the current dimension
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	void updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
};

}
}

#endif /* OPERATIONUPDOWNTESTLINEARBOUNDARY_HPP */
