/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "basis/linear/noboundary/algorithm_sweep/HierarchisationLinear.hpp"

namespace sg
{

namespace detail
{

HierarchisationLinear::HierarchisationLinear(GridStorage* storage) : storage(storage)
{
}

HierarchisationLinear::~HierarchisationLinear()
{
}

void HierarchisationLinear::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	rec(source, result, index, dim, 0.0, 0.0);
}

void HierarchisationLinear::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{
	// current position on the grid
	size_t seq = index.seq();
	// value in the middle, needed for recursive call and calculation of the hierarchical surplus
	double fm = source[seq];

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

	// hierarchisation
	result[seq] = fm - ((fl + fr)/2.0);
}

}	// namespace detail

}	// namespace sg
