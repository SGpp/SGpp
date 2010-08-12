/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONGAMMALINEAR_HPP
#define OPERATIONGAMMALINEAR_HPP

#include "algorithm/pde/UpDownTwoOpDims.hpp"

namespace sg
{

/**
 * Implements the Gamma-Operation (corresponds to matrix E in Master's thesis), that is needed
 * the solve the multidimensional Black Scholes
 * equation, on grids with fix Dirichlet-0-Boundaries.
 *
 * @version $HEAD$
 */
class OperationGammaLinear : public UpDownTwoOpDims
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param coef vector that contains the constant coefficients of this operation
	 */
	OperationGammaLinear(GridStorage* storage, DataMatrix& coef);

	/**
	 * Destructor
	 */
	virtual ~OperationGammaLinear();

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
	 * down-Gradient step in dimension <i>dim</i> applies the x phi dphi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downOpDimOne(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * up-Gradient step in dimension <i>dim</i> applies the x phi dphi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upOpDimOne(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * down-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downOpDimTwo(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * up-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upOpDimTwo(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * down-Gradient multiplied with a squared x step in dimension <i>dim</i> applies the x^2 dphi dphi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	void downOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim);

	/**
	 * up-Gradient multiplied with a squared x step in dimension <i>dim</i> applies the x^2 dphi dphi operation
	 * in one dimension
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	void upOpDimOneAndOpDimTwo(DataVector& alpha, DataVector& result, size_t dim);
};

}

#endif /* OPERATIONGAMMALINEAR_HPP */
