/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACEMODLINEAR_HPP
#define OPERATIONLAPLACEMODLINEAR_HPP

#include "algorithm/pde/UpDownOneOpDim.hpp"
using namespace sg::base;

namespace sg
{

/**
 * Implementation of Laplace for mod linear functions
 *
 * @version $HEAD$
 */
class OperationLaplaceModLinear : public UpDownOneOpDim
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationLaplaceModLinear(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~OperationLaplaceModLinear();

protected:
	virtual void up(DataVector& alpha, DataVector& result, size_t dim);

	virtual void down(DataVector& alpha, DataVector& result, size_t dim);

	virtual void downOpDim(DataVector& alpha, DataVector& result, size_t dim);

	virtual void upOpDim(DataVector& alpha, DataVector& result, size_t dim);
};

}

#endif /* OPERATIONLAPLACEMODLINEAR_HPP */
