/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include "basis/modlinear/algorithm_sweep/DehierarchisationModLinear.hpp"

namespace sg
{

namespace detail
{

DehierarchisationModLinear::DehierarchisationModLinear(GridStorage* storage) : storage(storage)
{
}

DehierarchisationModLinear::~DehierarchisationModLinear()
{
}

void DehierarchisationModLinear::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	rec(source, result, index, dim, 0.0, 0.0);
}

void DehierarchisationModLinear::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{
	// current position on the grid
	size_t seq = index.seq();
	// value in the middle, needed for recursive call and calculation of the hierarchical surplus
	double fm = source[seq];

	// dehierarchisation
	fm += ((fl + fr)/2.0);
	result[seq] = fm;

	GridStorage::index_type::level_type l;
	GridStorage::index_type::index_type i;

	index.get(dim, l, i);

	// recursive calls for the right and left side of the current node
	if(index.hint() == false)
	{
		double fltemp = fl;
		double frtemp = fr;

		// When we descend the hierarchical basis we have to modify the boundary values
		// in case the index is 1 or (2^l)-1 or we are on the first level
		// level 1, constant function
		if(l == 1)
		{
			// constant function
			fltemp = fm;
			frtemp = fm;
		}
		// left boundary
		else if(i == 1)
		{
			double ftemp;
			ftemp = fr - fm;
			fltemp = fm - ftemp;
		}
		// right boundary
		else if(static_cast<int>(i) == static_cast<int>((1 << l)-1))
		{
			double ftemp;
			ftemp = fl - fm;
			frtemp = fm - ftemp;
		}
		// inner functions
		else
		{
		}

		// descend left
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fltemp, fm);
		}

		// descend right
		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fm, frtemp);
		}

		// ascend
		index.up(dim);
	}
}

}	// namespace detail

}	// namespace sg
