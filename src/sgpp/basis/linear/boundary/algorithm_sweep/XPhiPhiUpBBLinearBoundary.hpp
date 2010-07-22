/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi (qic@in.tum.de) Stefanie Schraufstetter (schraufs@in.tum.de)

#ifndef PHIPHIUPBBLINEARBOUNDARY_HPP
#define PHIPHIUPBBLINEARBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"

namespace sg
{

namespace detail
{

/**
 * Implementation of sweep operator (): 1D Up for
 * Bilinearform \f$\int_{x} x \phi(x) \phi(x) dx\f$
 * on linear boundary grids
 */
class XPhiPhiUpBBLinearBoundary : public XPhiPhiUpBBLinear
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	XPhiPhiUpBBLinearBoundary(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~XPhiPhiUpBBLinearBoundary();

	/**
	 * This operations performs the calculation of up in the direction of dimension <i>dim</i>
	 *
	 * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
	 * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
	 * result)
	 * So please assure that both functions do exist!
	 *
	 * @param source DataVector that contains the gridpoint's coefficients (values from the vector of the laplace operation)
	 * @param result DataVector that contains the result of the up operation
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);
};

} // namespace detail

} // namespace sg

#endif /* PHIPHIUPBBLINEARBOUNDARY_HPP */
