/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONLAPLACELINEAR_HPP
#define OPERATIONLAPLACELINEAR_HPP

#include "algorithm/pde/UpDownOneOpDim.hpp"

namespace sg
{

/**
 * Implementation for linear functions of Laplace Operation, linear grids without boundaries
 *
 * @version $HEAD$
 */
class OperationLaplaceLinear: public UpDownOneOpDim
{
public:
	/**
	 * Construtor of OperationLaplaceLinear
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationLaplaceLinear(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~OperationLaplaceLinear();

protected:
#ifndef USEOMPTHREE
	virtual void specialOP(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim);
#endif

#ifdef USEOMPTHREE
	virtual void specialOP_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim);
#endif

	virtual void up(DataVector& alpha, DataVector& result, size_t dim);

	virtual void down(DataVector& alpha, DataVector& result, size_t dim);

	virtual void downOpDim(DataVector& alpha, DataVector& result, size_t dim);

	virtual void upOpDim(DataVector& alpha, DataVector& result, size_t dim);
};

}

#endif /* OPERATIONLAPLACELINEAR_HPP */
