/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONDELTALOGLINEAR_HPP
#define OPERATIONDELTALOGLINEAR_HPP

#include "algorithm/pde/UpDownOneOpDim.hpp"

namespace sg
{

/**
 * Implements the Delta-Operation, that is needed
 * the solve the multidimensional Black Scholes
 * equation
 *
 * This operation implements the Delta-Operation
 * for the log-transformed Black Scholes Equation
 * on grid's with 0-Dririchlet Boundaries.
 *
 * @version $HEAD$
 */
class OperationDeltaLogLinear: public UpDownOneOpDim
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param coef reference to a DataVector object that contains the bilinear form's constant coefficients
	 */
	OperationDeltaLogLinear(GridStorage* storage, DataVector& coef);

	/**
	 * Destructor
	 */
	virtual ~OperationDeltaLogLinear();

protected:
	/**
	 * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
	 * Applies the up-part of the one-dimensional mass matrix in one dimension.
	 * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
	 *
	 * @param dim dimension in which to apply the up-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void up(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
	 * Applies the down-part of the one-dimensional mass matrix in one dimension.
	 * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
	 *
	 * @param dim dimension in which to apply the down-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void down(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * down-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downOpDim(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * up-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upOpDim(DataVector& alpha, DataVector& result, size_t dim);
};

}

#endif /* OPERATIONDELTALOGLINEAR_HPP */
