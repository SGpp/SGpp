/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACELINEARBOUNDARY_HPP
#define OPERATIONLAPLACELINEARBOUNDARY_HPP

#include "algorithm/pde/UpDownOneOpDim.hpp"
using namespace sg::base;

namespace sg
{
namespace pde
{

/**
 * Implementation of Laplace for linear functions with boundaries
 *
 * @version $HEAD$
 */
class OperationLaplaceLinearBoundary: public UpDownOneOpDim
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationLaplaceLinearBoundary(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~OperationLaplaceLinearBoundary();


protected:
	virtual void up(DataVector& alpha, DataVector& result, size_t dim);

	virtual void down(DataVector& alpha, DataVector& result, size_t dim);

	virtual void downOpDim(DataVector& alpha, DataVector& result, size_t dim);

	virtual void upOpDim(DataVector& alpha, DataVector& result, size_t dim);
};

}
}

#endif /* OPERATIONLAPLACELINEARBOUNDARY_HPP */
