/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHCOARSENING_HPP
#define HASHCOARSENING_HPP

#include "data/DataVector.hpp"

#include "grid/GridStorage.hpp"
#include "grid/generation/CoarseningFunctor.hpp"

#include "exception/generation_exception.hpp"

#include <vector>
#include <cmath>
#include <utility>

#include <iostream>

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
	typedef std::pair<size_t, CoarseningFunctor::value_type> GridPoint;

	/**
	 * Performs coarsening on grid. It's possible to remove a certain number
	 * of gridpoints in one coarsening step. This number is specified within the
	 * declaration of the coarsening functor. Also the coarsening threshold is
	 * specified in the coarsening functor.
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a function used to determine if refinement is needed
	 * @param alpha pointer to the gridpoints' coefficients removed points must also be considered in this vector
	 */
	void free_coarsen(GridStorage* storage, CoarseningFunctor* functor, DataVector* alpha)
	{
		// chech if the grid has any points
		if(storage->size() == 0)
		{
			throw generation_exception("storage empty");
		}

		// Perepare temp-data in order to determine the removable grid points
		// -> leafs with minimal surplus
		size_t remove_num = functor->getRemovementsNum();

		// create an array that will contain the GridPoints
		// (pair of the grid Point's index and its surplus)
		// that should be removed
		GridPoint* removePoints = new GridPoint[remove_num];

		// init the removePoints array:
		// set initial surplus and set all indices to zero
		for (size_t i = 0; i<remove_num; i++){
			removePoints[i].second = functor->start();
			removePoints[i].first = 0;
		}

		// help variable to store the gridpoint with highest
		// surplus in removePoints
		size_t max_idx = 0;

		for(GridStorage::grid_map_iterator iter = storage->begin(); iter != storage->end(); iter++)
		{
			index_type* index = iter->first;

			if (index->isLeaf())
			{
				CoarseningFunctor::value_type current_value = (*functor)(storage, iter->second);
				if(current_value < removePoints[max_idx].second)
				{
					//Replace the maximum point array of removable candidates, find the new maximal point
					removePoints[max_idx].second = current_value;
					removePoints[max_idx].first = iter->second;

					// find new maximum entry
					max_idx = 0;
					for (size_t i = 1; i < remove_num; i++)
					{
						if(removePoints[i].second > removePoints[max_idx].second)
						{
							max_idx = i;
						}
					}
				}
			}
		}

		//DEBUG : print list of removable candidates
		//for (size_t i = 0; i < remove_num; i++)
		//{
		//	std::cout << "Index: " << removePoints[i].first << " with surplus " << removePoints[i].second << std::endl;
		//}
		//std::cout << std::endl;

		// remove the marked grid point if their surplus
		// is below the given threshold
		CoarseningFunctor::value_type threshold = functor->getCoarseningThreshold();
		CoarseningFunctor::value_type initValue = functor->start();

		// copy GridStorage in order to shrink grid
		GridStorage tempStorage(*storage);

		// vector to save remaining points
		std::vector<size_t> remainingIndex;

		// remove all grid points;
		storage->emptyStorage();

		// re-insert grid points
		for(size_t i = 0; i < tempStorage.size(); i++)
		{
			index_type* index = tempStorage[i];

			if (shouldRemovePoint(removePoints, i, initValue, threshold, remove_num) == false)
			{
				storage->insert(*index);
				remainingIndex.push_back(i);
			}
		}

		// DEBUG
		//std::cout << "List of remaining GridPoints (indices)" << std::endl;
		//for (size_t i = 0; i < remainingIndex.size(); i++)
		//{
		//	std::cout << remainingIndex[i] << " ";
		//}
		//std::cout << std::endl << std::endl;

		// make the grid a consistent grid again -> the leaf property was destroyed
		// by removing grid points -> re-calculate it.
		storage->recalcLeafProperty();

		// Drop Elements from DataVector
		alpha->restructure(remainingIndex);
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

		for(GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++)
		{
			index = *(iter->first);

			if (index.isLeaf())
			{
				counter++;
			}
		}

		return counter;
	}


protected:

	/**
	 * Determines if a gird point should be removed from the grid
	 *
	 * @param removePoints point to an array of pairs of index and surplus of candidates fro removing
	 * @param current_index current grid point which is tested
	 * @param initValue init value of the removable points array for the surpluses
	 * @param threshold threshold used in removing decision
	 * @param remove_num number of elements in removePoints
	 */
	bool shouldRemovePoint(GridPoint* removePoints, size_t current_index, CoarseningFunctor::value_type initValue, CoarseningFunctor::value_type threshold, size_t remove_num)
	{
		for (size_t i = 0; i < remove_num; i++)
		{
			if (removePoints[i].first == current_index)
			{
				if(removePoints[i].second < initValue && removePoints[i].second <= threshold)
				{
					return true;
				}
				else
				{
					return false;
				}
			}
		}

		return false;
	}
};

}

#endif /* HASHCOARSENING_HPP */
