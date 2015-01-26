/* ****************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
 *
 * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data points
 */
class OperationMultipleEval {
protected:
	Grid &grid;
	DataMatrix &dataset;
	bool isPrepared;

public:
	/**
	 * Constructor
	 *
	 * @param grid the sparse grid used for this operation
	 * @param dataset data set that should be evaluated on the sparse grid, a operation may create a copy of the dataset
	 */
	OperationMultipleEval(SGPP::base::Grid &grid, DataMatrix &dataset) :
			grid(grid), dataset(dataset), isPrepared(false) {
	}

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEval() {
	}

	/**
	 * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
	 *
	 * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void mult(DataVector& alpha, DataVector& result) = 0;

	/**
	 * Multiplication of @f$B@f$ with vector @f$\alpha@f$
	 *
	 * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual void multTranspose(DataVector& source, DataVector& result) = 0;

	/**
	 * Evaluate multiple datapoints with the specified grid
	 *
	 * @param alpha surplus vector of the grid
	 * @param result result of the evaluations
	 */
	void eval(DataVector &alpha, DataVector &result) {
		this->mult(alpha, result);
	}

	/**
	 * Used for kernel-specific setup like special data structures that are defined from the current state of
	 * the grid. This function is by default called with each "mult()", "multTranspose()" or evaluation operation
	 * and can be ignored from an external perspective.
	 * This is not overridden by every kernel.
	 */
	virtual void prepare() {
	}

	/**
	 * Tells the kernel that the grid hasn't changed since the last call to "prepare()". The kernel will therefore
	 * assume that additional "prepare()" calls can be skipped when "mult()" or "multTransposed()" or evaluation operations
	 * are performed.
	 * To be consistent, the call to "prepare()" has to be performed explicitly. Also, "setPrepared(false)" has to be called
	 * after a changed to the grid.
	 *
	 * @param isPrepared Tells the operation that is has been prepared for the current (or not)
	 */
	void setPrepared(bool isPrepared) {
		this->isPrepared = isPrepared;
	}

	/**
	 * Name of this implementation of the operation.
	 */
	virtual std::string getImplementationName() {
		throw new SGPP::base::operation_exception(
				"error: OperationMultipleEval::getImplementationName(): not implemented for this kernel");
	}
};

}
}
