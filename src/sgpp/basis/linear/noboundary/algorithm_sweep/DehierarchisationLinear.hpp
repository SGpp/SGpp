/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DEHIERARCHISATIONLINEAR_HPP
#define DEHIERARCHISATIONLINEAR_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

namespace detail
{

/**
 * Class that implements the dehierarchisation on a linear sparse grid. Therefore
 * the ()operator has to be implement in order to use the sweep algorithm for
 * the grid traversal
 */
class DehierarchisationLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// the grid object
	GridStorage* storage;

public:
	/**
	 * Constructor, must be bind to a grid
	 *
	 * @param storage the grid storage object of the the grid, on which the dehierarchisation should be executed
	 */
	DehierarchisationLinear(GridStorage* storage) : storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	virtual ~DehierarchisationLinear()
	{
	}

	/**
	 * Implements operator() needed by the sweep class during the grid traversal. This function
	 * is applied to the whole grid.
	 *
	 * @param source this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
	 * @param result this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	virtual void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		rec(source, result, index, dim, 0.0, 0.0);
	}

protected:

	/**
	 * Recursive dehierarchisaton algorithm, this algorithms works in-place -> source should be equal to result
	 *
	 * @todo add graphical explanation here
	 *
	 * @param source this DataVector holds the linear base coefficients of the sparse grid's ansatz-functions
	 * @param result this DataVector holds the node base coefficients of the function that should be applied to the sparse grid
	 * @param index a iterator object of the grid
	 * @param dim current fixed dimension of the 'execution direction'
	 * @param fl left value of the current region regarded in this step of the recursion
	 * @param fr right value of the current region regarded in this step of the recursion
	 */
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
	{
		// current position on the grid
		size_t seq = index.seq();
		// value in the middle, needed for recursive call and calculation of the hierarchical surplus
		double fm = source[seq];

		// dehierarchisation
		fm += ((fl + fr)/2.0);
		result[seq] = fm;

		// recursive calls for the right and left side of the current node
		if(index.hint() == false)
		{
			// descend left
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fl, fm);
			}

			// descend right
			index.step_right(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fm, fr);
			}

			// ascend
			index.up(dim);
		}
	}
};

}	// namespace detail

}	// namespace sg

#endif /* DEHIERARCHISATIONLINEAR_HPP */
