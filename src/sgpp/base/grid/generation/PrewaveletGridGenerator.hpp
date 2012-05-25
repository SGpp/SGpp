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

#ifndef PREWAVELETGRIDGENERATOR_HPP
#define PREWAVELETGRIDGENERATOR_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/GridGenerator.hpp"

namespace sg
{
namespace base
{

/**
 * GridGenerator for prewavelet grids without boundaries
 */
class PrewaveletGridGenerator: public GridGenerator
{
protected:
	/// pointer to the storage object
	GridStorage* storage;
	GridStorage* shadowstorage;
	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;

private:
	void addNeighbours(index_type& index, size_t current_dim,
			level_t target_level, GridStorage::grid_iterator& iter, GridStorage::grid_iterator& shadowIter);
	void insertParents(GridStorage::grid_iterator& iter, GridStorage::grid_iterator& shadowIter);
	void consolidateShadow();

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's storage object
	 */
	PrewaveletGridGenerator(GridStorage* storage, GridStorage* shadowstorage);

	/**
	 * Destructor
	 */
	virtual ~PrewaveletGridGenerator();

	virtual void regular(size_t level);
	virtual void full(size_t level);
	virtual void refine(RefinementFunctor* func);
	virtual size_t getNumberOfRefinablePoints();

	virtual void coarsen(CoarseningFunctor* func, DataVector* alpha);
	virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly);
	virtual size_t getNumberOfRemovablePoints();

	virtual void refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel);
	virtual size_t getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel);

};

}
}

#endif /* PREWAVELETGRIDGENERATOR_HPP */
