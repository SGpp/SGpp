/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef HASHREFINEMENTABSTRACT_HPP
#define HASHREFINEMENTABSTRACT_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/refinement_strategy/RefinementStrategy.hpp"

namespace sg
{
namespace base
{

//class RefinementStrategy;
/**
 * Standard free refinement class for sparse grids without boundaries
 */
class HashRefinementAbstract
{
public:
	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;

	/**
	 * Refines a grid according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if 
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	virtual void free_refine(GridStorage* storage, RefinementFunctor* functor)=0;


	/**
	 * Computes and returns the number of grid points, which can be refined. 
	 * This is the number of grid points that have at least one child missing.
	 *
	 * @param storage hashmap that stores the grid points
	 * @return The number of grid points that can be refined
	 */
	virtual size_t getNumberOfRefinablePoints(GridStorage* storage) = 0;

	virtual void strategy_refine(GridStorage* storage, RefinementStrategy& refinement_strategy);
	/*{
		refinement_strategy.refine(storage, this);
	}*/

	virtual void refine_gridpoint_1d(GridStorage * storage, index_type & index, size_t d) = 0;

	~HashRefinementAbstract(){};

protected:
	/**
	 * This method refines a grid point by generating the children in every dimension
	 * of the grid and all their missing ancestors by calling create_gridpoint().
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param refine_index The index in the hashmap of the point that should be refined
	 */
	virtual void refine_gridpoint(GridStorage* storage, size_t refine_index) = 0;

	/**
	 * This method creates a new point on the grid. It checks if some parents or
	 * children are needed in other dimensions.
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param index The point that should be inserted
	 */
	virtual void create_gridpoint(GridStorage* storage, index_type& index) = 0;


	/**
	 * Returns the index of the first occurrence of minimal element in array.
	 * Used to find which entry is to be replaced next searching the maximum ones.
	 *
	 * @param array array with values
	 * @param length length of array
	 *
	 * @return index of the first occurrence of minimal element in array
	 */
	virtual size_t getIndexOfMin(RefinementFunctor::value_type* array, size_t length);


	virtual void create_gridpoint_subroutine(GridStorage* storage, index_type& index){
		if(!storage->has_key(&index))
			{
				// save old leaf value
				bool saveLeaf = index.isLeaf();
				index.setLeaf(false);
				create_gridpoint(storage, index);
				// restore leaf value
				index.setLeaf(saveLeaf);
			}
			else
			{
				// set stored index to false
				(storage->get((storage->find(&index))->second))->setLeaf(false);
			}
	};

	virtual void create_gridpoint_1d(index_type& index,
			size_t d, GridStorage * storage, index_t& source_index, level_t& source_level);




};

}
}

#endif /* HASHREFINEMENTABSTRACT_HPP */
