/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#ifndef OPERATIONHESTONFLINEAR_HPP
#define OPERATIONHESTONFLINEAR_HPP

#include "pde/algorithm/UpDownOneOpDim.hpp"

namespace sg
{
namespace finance
{

/**
 * Todo: fix comments
 */
class OperationHestonFLinear: public sg::pde::UpDownOneOpDim
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's sg::base::GridStorage object
	 * @param coef reference to a sg::base::DataVector object that contains the bilinear form's constant coefficients
	 */
	OperationHestonFLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef);

	/**
	 * Destructor
	 */
	virtual ~OperationHestonFLinear();

protected:
	/**
	 * Todo: fix comments
	 */
	virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	/**
	 * Todo: fix comments
	 */
	virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	/**
	 * Todo: fix comments
	 */
	virtual void downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	/**
	 * Todo: fix comments
	 */
	virtual void upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
};

}
}

#endif /* OPERATIONHESTONFLINEAR_HPP */
