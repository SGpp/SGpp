/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#ifndef HASHREFINEMENTABSTRACT_HPP
#define HASHREFINEMENTABSTRACT_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/refinement_strategy/RefinementStrategy.hpp"

namespace sg
{
namespace base
{

/**
 * Abstract refinement class for sparse grids
 */
class AbstractRefinement
{
public:
    typedef GridStorage::index_type index_type;
    typedef GridStorage::index_pointer index_pointer;
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


    /**
     * Refine grid points using a specified strategy
     * @param storage hashmap that stores the grid points
     * @param index point to refine
     * @param d direction
     */
    virtual void strategy_refine(GridStorage* storage, RefinementStrategy& refinement_strategy);

    /**
     * Refine one grid point along a single direction
     * @param storage hashmap that stores the grid points
     * @param index point to refine
     * @param d direction
     */
    virtual void refine_gridpoint_1d(GridStorage * storage, index_type & index, size_t d) = 0;

    /**
     * Refine one grid point along a single direction
     * @param storage hashmap that stores the grid points
     * @param index point to refine
     * @param d direction
     */
    void refine_gridpoint_1d(GridStorage * storage, index_pointer index, size_t d) {
      refine_gridpoint_1d(storage, *index, d);
    }

    bool is_refinable(GridStorage* storage, index_type& index);

    /**
     * Destructor
     */
    ~AbstractRefinement()
    {}
    ;

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
     * Subroutine for grid point creation.
     *
     * @param storage hashmap that stores the gridpoints
     * @param index The point that should be inserted
     */
    virtual void create_gridpoint_subroutine(GridStorage* storage, index_type& index)
    {
        // For efficiency this function is defined the header file, this way it
        // be easily inlined by compiler.
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


    /**
     * Creates children grid points along single direction
     *
     * @param index The point that should be refined
     * @param d direction
     * @param storage hashmap that stores the gridpoints
     * @param source_index index value in the dimension d
     * @param source_level level value in the dimension d
     */
    virtual void create_gridpoint_1d(index_type& index,
                                     size_t d, GridStorage * storage,
                                     index_t& source_index, level_t& source_level);




};

}
}

#endif /* HASHREFINEMENTABSTRACT_HPP */
