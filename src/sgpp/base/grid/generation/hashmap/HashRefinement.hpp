/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef HASHREFINEMENT_HPP
#define HASHREFINEMENT_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"

namespace sg
{
namespace base
{

/**
 * Abstract free refinement class for sparse grids
 */
class HashRefinement: public HashRefinementAbstract
{

public:


	/**
	 * Refines a grid according to a RefinementFunctor provided.
	 * Refines up to RefinementFunctor::getRefinementsNum() grid points if 
	 * possible, and if their refinement value is larger than RefinementFunctor::start()
	 * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a RefinementFunctor specifying the refinement criteria
	 */
	void free_refine(GridStorage* storage, RefinementFunctor* functor);


	/**
	 * Computes and returns the number of grid points, which can be refined. 
	 * This is the number of grid points that have at least one child missing.
	 *
	 * @param storage hashmap that stores the grid points
	 * @return The number of grid points that can be refined
	 */
	size_t getNumberOfRefinablePoints(GridStorage* storage);


	void refine_gridpoint_1d(GridStorage * storage, index_type & index, size_t d);

protected:
	/**
	 * This method refines a grid point by generating the children in every dimension
	 * of the grid and all their missing ancestors by calling create_gridpoint().
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param refine_index The index in the hashmap of the point that should be refined
	 */
	void refine_gridpoint(GridStorage* storage, size_t refine_index);
	/**
	 * This method creates a new point on the grid. It checks if some parents or
	 * children are needed in other dimensions.
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param index The point that should be inserted
	 */
	void create_gridpoint(GridStorage* storage, index_type& index);

	/*void refine_gridpoint_1d(index_type& index,
				size_t d, GridStorage * storage, index_t& souce_index, level_t& source_level);*/


};

}
}

#endif /* HASHREFINEMENT_HPP */
