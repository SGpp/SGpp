/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DEHIERARCHISATIONLINEARBOUNDARY_HPP
#define DEHIERARCHISATIONLINEARBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "basis/linear/noboundary/algorithm_sweep/DehierarchisationLinear.hpp"

namespace sg
{

namespace detail
{

/**
 * Class that implements the dehierarchisation on a linear sparse grid with boundaries. Therefore
 * the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal
 */
class DehierarchisationLinearBoundary : public DehierarchisationLinear
{
public:
	/**
	 * Constructor, must be bind to a grid
	 *
	 * @param storage the grid storage object of the the grid, on which the dehierarchisation should be executed
	 */
	DehierarchisationLinearBoundary(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~DehierarchisationLinearBoundary();

	/**
	 * Implements operator() needed by the sweep class during the grid traversal. This function
	 * is applied to the whole grid.
	 *
	 * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
	 * If one is missing this code might produce some bad errors (segmentation fault, wrong calculation
	 * result)
	 * So please assure that both functions do exist!
	 *
	 * @param source this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
	 * @param result this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim);
};

}	// namespace detail

}	// namespace sg

#endif /* DEHIERARCHISATIONLINEARBOUNDARY_HPP */
