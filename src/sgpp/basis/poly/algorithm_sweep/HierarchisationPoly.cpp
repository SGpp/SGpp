/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Kilian Roehner (roehner@tum.de)

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include <cmath>

#include "basis/poly/algorithm_sweep/HierarchisationPoly.hpp"

namespace sg
{

namespace base
{

HierarchisationPoly::HierarchisationPoly(GridStorage* storage, SPolyBase base) : storage(storage), base(base)
{
}

HierarchisationPoly::~HierarchisationPoly()
{
}

void HierarchisationPoly::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	// double &koeffs;
	DataVector koeffs(index.getGridDepth(dim));
	koeffs.setAll(0.0);
	rec(source, result, index, dim, koeffs);
}

void HierarchisationPoly::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, DataVector& koeffs)
{
	// current position on the grid
	size_t seq = index.seq();
	
	// value in the middle, needed for recursive call and calculation of the hierarchical surplus
	double fm = source[seq];
	
	GridStorage::index_type::level_type cur_lev;
	GridStorage::index_type::index_type cur_ind;
	
	// get current level and index from grid
	index.get(seq, cur_lev, cur_ind);
	
	// calculate the current absolute position
	double abs = cur_ind*(pow(2.0, -1*cur_lev));
	
	// calculate the weighted "so far" sum of the basis funktions at this point
	double val = this->base.evalHierToTop(cur_lev, cur_ind, abs, koeffs);
	
	// hierarchisation
	result[seq] = fm - val;
	
	// recursive calls for the right and left side of the current node
	if(index.hint() == false)
	{
		koeffs[cur_lev] = result[seq];
	
		// descend left
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, koeffs);
		}

		// descend right
		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, koeffs);
		}

		// ascend
		index.up(dim);
		
		koeffs[cur_lev] = 0.0;
	}
	
}

}	// namespace base

}	// namespace sg
