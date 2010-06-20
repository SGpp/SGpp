/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHCOARSENING_HPP
#define HASHCOARSENING_HPP

#include "grid/GridStorage.hpp"
#include "grid/generation/CoarseningFunctor.hpp"

#include "exception/generation_exception.hpp"

#include <vector>
#include <cmath>

namespace sg
{

/**
 * Standard free coarsening class for sparse grids without boundaries
 */
class HashCoarsening
{
public:
	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;

	/**
	 * Performs coarsening on grid. It's possible to remove a certain number
	 * of gridpoints in one coarsening step. This number is specified within the
	 * declaration of the coarsening functor. Also the coarsening threshold is
	 * specified in the coarsening functor.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a function used to determine if refinement is needed
	 */
	void free_coarsen(GridStorage* storage, CoarseningFunctor* functor)
	{
		// chech if the grid has any points
		if(storage->size() == 0)
		{
			throw generation_exception("storage empty");
		}

		// Perepare temp-data in order to determine the removable grid points
		// -> leafs with minimal surplus
		size_t remove_num = functor->getRemovementsNum();
		RefinementFunctor::value_type* min_values = new RefinementFunctor::value_type[remove_num];
		size_t* min_indexes = new size_t [refinements_num];

		// init temp data
		for (size_t i = 0; i<remove_num; i++){
			min_values[i] = functor->start();
			min_indexes[i] = 0;
		}
		size_t max_idx = 0;

		RefinementFunctor::value_type min_value = min_values[max_idx];
		size_t min_index = min_indexes[max_idx];

		index_type index;
		GridStorage::grid_map_iterator end_iter = storage->end();

		for(GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++)
		{
			index = *(iter->first);

			GridStorage::grid_map_iterator child_iter;

			if (index->isLeaf())
			{
				RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
				if(current_value < min_value)
				{
					//Replace the maximum point array of removable candidates, find the new maximal point
					min_values[max_idx] = current_value;
					min_indexes[max_idx] = iter->second;
					max_idx = getIndexOfMax(min_values, remove_num);
					min_value = min_values[max_idx];
					break;
				}
			}
		}


		//can refine grid on several points
		double threshold = functor->getCoarseningThreshold();
		for (size_t i = 0; i < refinements_num; i++){
			max_value = max_values[i];
			max_index = max_indexes[i];

			if(max_value > functor->start() && fabs(max_value) >= threshold)
			{
				remove_gridpoint(storage, max_index);
			}
		}

		delete[] max_values;
		delete[] max_indexes;

		// make the grid a consistent grid again -> the leaf property was destroyed
		// by removing grid points -> re-calculate it.
		storage->recalcLeafProperty();
	}

	/**
	 * Calculates the number of points, which can be refined
	 *
	 * @param storage hashmap that stores the grid points
	 */
	size_t getNumberOfRemovablePoints(GridStorage* storage)
	{
		size_t counter = 0;

		if(storage->size() == 0)
		{
			throw generation_exception("storage empty");
		}

		index_type index;
		GridStorage::grid_map_iterator end_iter = storage->end();

		// I think this may be depedent on local support
		for(GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++)
		{
			index = *(iter->first);

			if (index->isLeaf())
			{
				counter++;
			}
		}

		return counter;
	}


protected:
	/**
	 * This method refines a grid point be generating the children in every dimension
	 * of the grid.
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param refine_index the index in the hashmap of the point that should be refined
	 */
	void remove_gridpoint(GridStorage* storage, size_t refine_index)
	{
		index_type index((*storage)[refine_index]);

		//Sets leaf property of index, which is refined to false
		//(*storage)[refine_index]->setLeaf(false);

	}

	/**
	 * Returns the index of the first accurance of minimal element in array
	 *
	 * @param array array with ???
	 * @param length length of ??
	 *
	 * @return index of the first accurance of minimal element in array
	 */
	size_t getIndexOfMax(RefinementFunctor::value_type* array, size_t length)
	{
		size_t max_idx = 0;
		for (size_t i = 1; i < length; i++)
		{
			if(array[i] > array[max_idx])
				max_idx = i;
		}

		return max_idx;
	}
};

}

#endif /* HASHCOARSENING_HPP */
