/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#ifndef OPERATIONHESTONBLINEAR_HPP
#define OPERATIONHESTONBLINEAR_HPP

#include "pde/algorithm/UpDownOneOpDim.hpp"

namespace sg
{
namespace finance
{

/**
 * Implements the Heston B-Operation (corresponds to matrix B in Master's thesis), that is needed
 * the solve the multidimensional Heston
 * equation, on grids with fix Dirichlet-0-Boundaries.
 *
 * @version $HEAD$
 */
class OperationHestonBLinear : public sg::pde::UpDownOneOpDim
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's sg::base::GridStorage object
	 * @param coef reference to a sg::base::DataVector object that contains the bilinear form's constant coefficients
	 */
	OperationHestonBLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef);

	/**
	 * Destructor
	 */
	virtual ~OperationHestonBLinear();

protected:
	/**
	 * Todo: improve comments to match existing format
	 *
	 * @param dim dimension in which to apply the up-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	/**
	 * Todo: improve comments to match existing format
	 *
	 * @param dim dimension in which to apply the down-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	/**
	 * Todo: improve comments to match existing format
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	/**
	 * Todo: improve comments to match existing format
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
};

}
}

#endif /* OPERATIONDELTALINEAR_HPP */
