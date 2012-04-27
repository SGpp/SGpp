/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHREFINEMENTBOUNDARIES_HPP
#define HASHREFINEMENTBOUNDARIES_HPP

#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"


namespace sg
{
namespace base
{

/**
 * Standard free refinement class for sparse grids with boundaries
 */
class HashRefinementBoundaries: public HashRefinementAbstract
{
public:

	/**
	 * Performs the refinement on grid
	 *
	 * @param storage hashmap that stores the grid points
	 * @param functor a function used to determine if refinement is needed
	 */
	void free_refine(GridStorage* storage, RefinementFunctor* functor);


	/**
	 * Calculates the number of points, which can be refined
	 *
	 * @param storage hashmap that stores the grid points
	 */
	size_t getNumberOfRefinablePoints(GridStorage* storage);

	void refine_gridpoint_1d(GridStorage * storage, HashRefinementAbstract::index_type & index, size_t d);

protected:
	/**
	 * This method refines a grid point be generating the children in every dimension
	 * of the grid.
	 *
	 * @param storage hashmap that stores the gridpoints
	 * @param refine_index the index in the hashmap of the point that should be refined
	 */
	void refine_gridpoint(GridStorage* storage, size_t refine_index);


	/**
	 * Wrapper for the two functions create_gridpoint_general and
	 * create_gridpoint_levelZeroConsistency which have both to be
	 * executed if a gridpoint is refined
	 *
	 * @param storage hashmap that stores the gridpoinrs
	 * @param index the point that should be inserted
	 */
	void create_gridpoint(GridStorage* storage, index_type& index);


	/**
	 * This method creates a new point on the grid. It checks if some parents or
	 * children are needed in other dimensions.
	 *
	 * @param storage hashmap that stores the gridpoinrs
	 * @param index the point that should be inserted
	 */
	void create_gridpoint_general(GridStorage* storage, index_type& index);


	/**
	 * Assures that we have always a consistent grid with both functions
	 * 0,0 and 0,1 on level zero
	 *
	 * @param storage hashmap that stores the gridpoinrs
	 * @param index the point that should be inserted
	 */
	void create_gridpoint_levelZeroConsistency(GridStorage* storage, index_type& index);

	/**
         * Creates children grid points along single direction
         *
         * @param index The point that should be refined
         * @param d direction
         * @param storage hashmap that stores the gridpoints
         * @param source_index index value in the dimension d
         * @param source_level level value in the dimension d
         */
	void create_gridpoint_1d(index_type& index,
				 size_t d, GridStorage * storage,
				 index_t& souce_index, level_t& source_level);




};

}
}

#endif /* HASHREFINEMENTBOUNDARIES_HPP */
