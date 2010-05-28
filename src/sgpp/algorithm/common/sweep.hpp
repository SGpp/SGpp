/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SWEEP_HPP
#define SWEEP_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include <vector>
#include <utility>
#include <iostream>

namespace sg
{

/**
 * Standard sweep operation
 * FUNC should be a class with overwritten operator(). For an example see laplace_up_functor in laplace.hpp.
 * It must be default constructable or copyable.
 * STORAGE must provide a grid_iterator supporting left_child, step_right, up, hint and seq.
 */
template<class FUNC>
class sweep
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Object of FUNC, this is executed by sweep
	FUNC functor;
	/// Pointer to the grid's storage object
	GridStorage* storage;

public:
	/**
	 * Create a new sweep object with a default constructed functor
	 *
	 * @param storage the storage that contains the grid points
	 */
	sweep(GridStorage* storage) : functor(), storage(storage)
	{
	}

	/**
	 * Create a new sweep object with a copied functor
	 *
	 * @param functor the functor that is executed on the grid
	 * @param storage the storage that contains the grid points
	 */
	sweep(FUNC& functor, GridStorage* storage) : functor(functor), storage(storage)
	{
	}

	/**
	 * Destructor
	 */
	~sweep()
	{
	}

	/**
	 * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep
	 * Boundaries are not regarded
	 *
	 * @param source a DataVector containing the source coefficients of the grid points
	 * @param result a DataVector containing the result coefficients of the grid points
	 * @param dim_sweep the dimension in which the functor is executed
	 */
	void sweep1D(DataVector& source, DataVector& result, size_t dim_sweep)
	{
		// generate a list of all dimension (-dim_sweep) from dimension recursion unrolling
		std::vector<size_t> dim_list;
		for(size_t i = 0; i < storage->dim(); i++)
		{
			if(i != dim_sweep)
			{
				dim_list.push_back(i);
			}
		}

		grid_iterator index(storage);

		sweep_rec(source, result, index, dim_list, storage->dim()-1, dim_sweep);
	}

	/**
	 * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep
	 * Boundaries are regarded
	 *
	 * @param source a DataVector containing the source coefficients of the grid points
	 * @param result a DataVector containing the result coefficients of the grid points
	 * @param dim_sweep the dimension in which the functor is executed
	 */
	void sweep1D_Boundary(DataVector& source, DataVector& result, size_t dim_sweep)
	{
		// generate a list of all dimension (-dim_sweep) from dimension recursion unrolling
		std::vector<size_t> dim_list;
		for(size_t i = 0; i < storage->dim(); i++)
		{
			if(i != dim_sweep)
			{
				dim_list.push_back(i);
			}
		}

		grid_iterator index(storage);
		index.resetToLevelZero();

		sweep_Boundary_rec(source, result, index, dim_list, storage->dim()-1, dim_sweep);
	}

protected:

	/**
	 * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep.
	 * Boundaries are not regarded
	 *
	 * @param source coefficients of the sparse grid
	 * @param result coefficients of the function computed by sweep
	 * @param index current grid position
	 * @param dim_list list of dimensions, that should be handled
	 * @param dim_rem number of remaining dims
	 * @param dim_sweep static dimension, in this dimension the functor is executed
	 */
	void sweep_rec(DataVector& source, DataVector& result, grid_iterator& index,
				std::vector<size_t>& dim_list, size_t dim_rem, size_t dim_sweep)
	{
		functor(source, result, index, dim_sweep);

		// dimension recursion unrolled
		for(size_t d = 0; d < dim_rem; d++)
		{
			size_t current_dim = dim_list[d];

			if(index.hint())
			{
				continue;
			}

			index.left_child(current_dim);
			if(!storage->end(index.seq()))
			{
				sweep_rec(source, result, index, dim_list, d+1, dim_sweep);
			}

			index.step_right(current_dim);
			if(!storage->end(index.seq()))
			{
				sweep_rec(source, result, index, dim_list, d+1, dim_sweep);
			}

			index.up(current_dim);
		}
	}

	/**
	 * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep.
	 * Boundaries are regarded
	 *
	 * @param source coefficients of the sparse grid
	 * @param result coefficients of the function computed by sweep
	 * @param index current grid position
	 * @param dim_list list of dimensions, that should be handled
	 * @param dim_rem number of remaining dims
	 * @param dim_sweep static dimension, in this dimension the functor is executed
	 */
	void sweep_Boundary_rec(DataVector& source, DataVector& result, grid_iterator& index,
				std::vector<size_t>& dim_list, size_t dim_rem, size_t dim_sweep)
	{
		if (dim_rem == 0)
		{
			functor(source, result, index, dim_sweep);
		}
		else
		{
			typedef GridStorage::index_type::level_type level_type;
			typedef GridStorage::index_type::index_type index_type;

			level_type current_level;
			index_type current_index;

			index.get(dim_list[dim_rem-1], current_level, current_index);

			// handle level greater zero
			if (current_level > 0)
			{
				// given current point to next dim
				sweep_Boundary_rec(source, result, index, dim_list, dim_rem-1, dim_sweep);

				if (!index.hint())
				{
					index.left_child(dim_list[dim_rem-1]);
					if(!storage->end(index.seq()))
					{
						sweep_Boundary_rec(source, result, index, dim_list, dim_rem, dim_sweep);
					}

					index.step_right(dim_list[dim_rem-1]);
					if(!storage->end(index.seq()))
					{
						sweep_Boundary_rec(source, result, index, dim_list, dim_rem, dim_sweep);
					}

					index.up(dim_list[dim_rem-1]);
				}
			}
			// handle level zero
			else
			{
				sweep_Boundary_rec(source, result, index, dim_list, dim_rem-1, dim_sweep);

				index.right_levelzero(dim_list[dim_rem-1]);
				sweep_Boundary_rec(source, result, index, dim_list, dim_rem-1, dim_sweep);

				if (!index.hint())
				{
					index.top(dim_list[dim_rem-1]);
					if(!storage->end(index.seq()))
					{
						sweep_Boundary_rec(source, result, index, dim_list, dim_rem, dim_sweep);
					}
				}

				index.left_levelzero(dim_list[dim_rem-1]);
			}
		}
	}
};

}

#endif /* SWEEP_HPP */
