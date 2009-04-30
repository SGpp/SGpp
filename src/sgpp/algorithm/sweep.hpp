/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef SWEEP_HPP
#define SWEEP_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

#include <vector>
#include <utility>
#include <iostream>

#ifdef USEOMPTEST
#include <omp.h>
#endif

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
#ifdef USEOMPTEST
	/// Number of threads max. used by sweep
	size_t max_threads;
	/// Counter that counts to current running threads
	size_t run_threads;
	/// Mutex to control the access on run_threads
	omp_lock_t run_lock;
#endif

public:
	/**
	 * Create a new sweep object with a default constructed functor
	 *
	 * @param storage the storage that contains the grid points
	 */
	sweep(GridStorage* storage) : functor(), storage(storage)
	{
#ifdef USEOMPTEST
		max_threads = 8;
		run_threads = 0;
		omp_init_lock(&run_lock);
#endif
	}

	/**
	 * Create a new sweep object with a copied functor
	 *
	 * @param functor the functor that is executed on the grid
	 * @param storage the storage that contains the grid points
	 */
	sweep(FUNC& functor, GridStorage* storage) : functor(functor), storage(storage)
	{
#ifdef USEOMPTEST
		max_threads = 8;
		run_threads = 0
		omp_init_lock(&run_lock);
#endif
	}

	/**
	 * Destructor
	 */
	~sweep()
	{
#ifdef USEOMPTEST
		omp_destroy_lock(&run_lock);
#endif
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

#ifdef USEOMPTEST
		max_threads = 8;
		//omp_set_nested(1);
		sweep_rec_parallel(source, result, index, dim_list, storage->dim()-1, dim_sweep);
		//omp_set_nested(0);
#else
		sweep_rec(source, result, index, dim_list, storage->dim()-1, dim_sweep);
#endif
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
	 * @todo check if it's possible to write a parallel implementation using OMP 3
	 *
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

#ifdef USEOMPTEST
	/**
	 * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep.
	 * Boundaries are regarded
	 *
	 * This Version is parallelized using the OpenMP 2 Sections part 1
	 *
	 * @param source coefficients of the sparse grid
	 * @param result coefficients of the function computed by sweep
	 * @param index current grid position
	 * @param dim_list list of dimensions, that should be handled
	 * @param dim_rem number of remaining dims
	 * @param dim_sweep static dimension, in this dimension the functor is executed
	 */
	void sweep_rec_parallel(DataVector& source, DataVector& result, grid_iterator& index,
				std::vector<size_t>& dim_list, size_t dim_rem, size_t dim_sweep)
	{
		if (dim_rem == 0)
		{
			functor(source, result, index, dim_sweep);
		}
		else
		{
			if (!index.hint())
			{
				if (run_threads < max_threads)
				{
					grid_iterator indexLeft(index);
					grid_iterator indexRight(index);
					increaseRunThreads(2);

					#pragma omp parallel sections num_threads(2) firstprivate(dim_rem, dim_sweep)
					{
						#pragma omp section
						{
							indexLeft.left_child(dim_list[dim_rem-1]);
							if(!storage->end(indexLeft.seq()))
							{
								sweep_rec_parallel(source, result, indexLeft, dim_list, dim_rem, dim_sweep);
							}
						}

						#pragma omp section
						{
							indexRight.right_child(dim_list[dim_rem-1]);
							if(!storage->end(indexRight.seq()))
							{
								sweep_rec_parallel(source, result, indexRight, dim_list, dim_rem, dim_sweep);
							}
						}
					}

					decreaseRunThreads(2);
				}
				else
				{
					index.left_child(dim_list[dim_rem-1]);
					if(!storage->end(index.seq()))
					{
						sweep_rec_parallel(source, result, index, dim_list, dim_rem, dim_sweep);
					}

					index.step_right(dim_list[dim_rem-1]);
					if(!storage->end(index.seq()))
					{
						sweep_rec_parallel(source, result, index, dim_list, dim_rem, dim_sweep);
					}

					index.up(dim_list[dim_rem-1]);
				}
			}

			// given current point to next dim
			sweep_rec_parallel(source, result, index, dim_list, dim_rem-1, dim_sweep);
		}
	}
#endif

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

#ifdef USEOMPTEST
	/**
	 * Increases the number of running threads for the current sweep instance
	 *
	 * \param add the delta that is added to run_threads
	 */
	void increaseRunThreads(size_t add)
	{
		omp_set_lock(&run_lock);
		run_threads += add;
		omp_unset_lock(&run_lock);
	}

	/**
	 * Decreases the number of running threads for the current sweep instance
	 *
	 * \param sub the delta that is subtracted from run_threads
	 */
	void decreaseRunThreads(size_t sub)
	{
		omp_set_lock(&run_lock);
		run_threads -= sub;
		omp_unset_lock(&run_lock);
	}
#endif
};

}

#endif /* SWEEP_HPP */
