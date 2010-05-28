/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/linear/noboundary/algorithm_sweep/DehierarchisationLinear.hpp"

namespace sg
{

namespace detail
{

DehierarchisationLinear::DehierarchisationLinear(GridStorage* storage) : storage(storage)
{
}

DehierarchisationLinear::~DehierarchisationLinear()
{
}

void DehierarchisationLinear::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	rec(source, result, index, dim, 0.0, 0.0);
}

void DehierarchisationLinear::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
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

}	// namespace detail

}	// namespace sg
