/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2010 Richard Röttger                                        */
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

#include "base/grid/generation/PrewaveletGridGenerator.hpp"
#include "base/grid/GridStorage.hpp"

#include "base/exception/generation_exception.hpp"

#include "base/grid/generation/hashmap/HashCoarsening.hpp"
#include "base/grid/generation/hashmap/HashRefinement.hpp"
#include "base/grid/generation/hashmap/HashGenerator.hpp"

#include <iostream>
//
namespace sg
{
namespace base
{

/**
 * An adaptive grid with prewavelet ansatz functions requires for operations
 * using the up-down algorithm shadow points. These shadow points a needed just
 * for data transport, thus they do not have an influence on the final function.
 * Please refer to sg::pde::UpDownOneOpDimWithShadow for more information.
 */
PrewaveletGridGenerator::PrewaveletGridGenerator(GridStorage* storage,
		GridStorage* shadowstorage) :
	storage(storage), shadowstorage(shadowstorage)
{
}

PrewaveletGridGenerator::~PrewaveletGridGenerator()
{
}

void PrewaveletGridGenerator::regular(size_t level)
{
	HashGenerator gen;
	gen.regular(this->storage, level);
}

void PrewaveletGridGenerator::full(size_t level)
{
	HashGenerator gen;
	gen.full(this->storage, level);
}

/**
 * Refines the grid and updates the shadow storage.
 */
void PrewaveletGridGenerator::refine(RefinementFunctor* func)
{
	HashRefinement refine;
	size_t start = this->storage->size();
	refine.free_refine(this->storage, func);
	size_t end = this->storage->size();
	//All added gridpoint are between [start,end[

	//Check if a gridpoint within the shadow storage is now part of the actual grid!
	for (size_t i = start; i < end; i++)
	{
		if (shadowstorage->find(storage->get(i)) != shadowstorage->end())
		{
			consolidateShadow();
			break;
		}
	}

	//Now add all missing neigbours to the shadowStorage
	for (size_t i = start; i < end; i++)
	{
		GridStorage::index_pointer index = this->storage->get(i);

		level_t sum = 0;
		for (size_t d = 0; d < storage->dim(); ++d)
		{
			index_t current_index;
			level_t current_level;
			index->get(d, current_level, current_index);
			sum += current_level;
		}

		GridStorage::grid_iterator iter(storage);
		GridStorage::grid_iterator shadowIter(shadowstorage);
		addNeighbours(*this->storage->get(i), 0, sum, iter, shadowIter);
	}

}

size_t PrewaveletGridGenerator::getNumberOfRefinablePoints()
{
	HashRefinement refine;
	return refine.getNumberOfRefinablePoints(this->storage);
}

/**
 * This function ensures that the special adaptive prewavelet grid points have parents.
 */
void PrewaveletGridGenerator::insertParents(GridStorage::grid_iterator& iter,
		GridStorage::grid_iterator& shadowIter)
{
	// Call parents in every dimension
	for (size_t d = 0; d < storage->dim(); d++)
	{
		index_t current_index;
		level_t current_level;
		iter.get(d, current_level, current_index);
		if (current_level == 1)
			continue;
		iter.up(d);
		shadowIter.up(d);
		this->storage->get(iter.seq())->setLeaf(false);

		//Ok, point is neither in storage, nor in shadowstorage ...
		if (storage->end(iter.seq()) && shadowstorage->end(shadowIter.seq()))
		{
			GridStorage::index_pointer new_index = new GridStorage::index_type(
					storage->dim());

			for (size_t dim = 0; dim < storage->dim(); ++dim)
			{
				index_t target_index;
				level_t target_level;

				iter.get(dim, target_level, target_index);

				new_index->set(dim, target_level, target_index);

			}
			shadowstorage->insert(*new_index);
			delete new_index;
			insertParents(iter, shadowIter);
		}
		iter.set(d, current_level, current_index);
		shadowIter.set(d, current_level, current_index);
	}
}

/**
 * For the shadow storage, the two left and two right neighbors in each
 * dimension of the refined point are required. This function only adds
 * the point which are not in the actual grid to the shadow storage.
 *
 * @param index point added during refinement
 */
void PrewaveletGridGenerator::addNeighbours(index_type& index,
		size_t current_dim, level_t target_level,
		GridStorage::grid_iterator& iter,
		GridStorage::grid_iterator& shadowIter)
{

	level_t sum = 0;
	for (size_t d = 0; d < storage->dim(); ++d)
	{
		index_t current_index;
		level_t current_level;
		iter.get(d, current_level, current_index);
		sum += current_level;
	}

	if (sum == target_level)
	{
		GridStorage::index_pointer new_index = new GridStorage::index_type(
				storage->dim());

		if (storage->end(iter.seq()) && shadowstorage->end(shadowIter.seq()))
		{
			//Ok, point is neither in storage, nor in shadowstorage ...
			//check if the border of index and iter touching
			for (size_t d = 0; d < storage->dim(); ++d)
			{
				index_t target_index;
				level_t target_level;

				index_t current_index;
				level_t current_level;

				iter.get(d, current_level, current_index);
				index.get(d, target_level, target_index);

				new_index->set(d, current_level, current_index);

				// The index cast to int is required to allow a negative index
				int target_left = (1.0 / (1 << target_level))
						* static_cast<double> (target_index - 3);
				int target_right = (1.0 / (1 << target_level))
						* static_cast<double> (target_index + 3);
				int current_left = (1.0 / (1 << current_index))
						* static_cast<double> (current_level + 3);
				int current_right = (1.0 / (1 << current_index))
						* static_cast<double> (current_level + 3);

				if (!(current_right > target_left || current_left
						< target_right))
				{
					return;
					delete new_index;
				}

			}

			//Yepp, the supports touching each other! Add point to shadow!
			shadowstorage->insert(*new_index);
			delete new_index;
			//Call for parents
			insertParents(iter, shadowIter);
		}

		return;
	}
	else if (sum > target_level)
	{
		return;
	}

	for (size_t d = current_dim; d < storage->dim(); d++)
	{
		index_t save_index;
		level_t save_level;
		iter.get(d, save_level, save_index); //Save current index

		iter.left_child(d);
		shadowIter.left_child(d);
		addNeighbours(index, d, target_level, iter, shadowIter);
		iter.set(d, save_level, save_index); //reset index
		shadowIter.set(d, save_level, save_index); //reset index

		iter.right_child(d);
		shadowIter.right_child(d);
		addNeighbours(index, d, target_level, iter, shadowIter);
		iter.set(d, save_level, save_index); //reset index
		shadowIter.set(d, save_level, save_index);
	}

}

/**
 * If during the refinement one or more points of the shadow register are added
 * to the actual grid then we have to remove these points from the shadow storage.
 */
void PrewaveletGridGenerator::consolidateShadow()
{

	GridStorage* temp = new GridStorage(storage->dim());
	for (size_t i = 0; i < shadowstorage->size(); i++)
	{
		temp->insert(*shadowstorage->get(i));
	}
	shadowstorage->emptyStorage();

	for (size_t i = 0; i < temp->size(); i++)
	{
		if (storage->find(temp->get(i)) == storage->end())
		{
			shadowstorage->insert(*temp->get(i));
		}
	}
	delete temp;
}

void PrewaveletGridGenerator::coarsen(CoarseningFunctor* func,
		DataVector* alpha)
{
	HashCoarsening coarsen;
	coarsen.free_coarsen(this->storage, func, alpha);
}

void PrewaveletGridGenerator::coarsenNFirstOnly(CoarseningFunctor* func,
		DataVector* alpha, size_t numFirstOnly)
{
	HashCoarsening coarsen;
	coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
}

size_t PrewaveletGridGenerator::getNumberOfRemovablePoints()
{
	HashCoarsening coarsen;
	return coarsen.getNumberOfRemovablePoints(this->storage);
}

void PrewaveletGridGenerator::refineMaxLevel(RefinementFunctor* func,
		unsigned int maxLevel)
{
}

size_t PrewaveletGridGenerator::getNumberOfRefinablePointsToMaxLevel(
		unsigned int maxLevel)
{
	return 0;
}
}
}
