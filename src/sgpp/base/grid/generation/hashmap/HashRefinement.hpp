/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef HASHREFINEMENT_HPP
#define HASHREFINEMENT_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"

namespace sg {
  namespace base {

    /**
     * Free refinement class for sparse grids
     */
    class HashRefinement: public AbstractRefinement {

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

        /**
         * Refine one grid point along a single direction
         * @param storage hashmap that stores the grid points
         * @param index point to refine
         * @param d direction
         */
        void refineGridpoint1D(GridStorage* storage, index_type& index, size_t d);

        void refineGridpoint1D(GridStorage* storage, HashGridIndex< unsigned int, unsigned int >* index, size_t d) {
          refineGridpoint1D(storage, *index, d);
        }

      protected:
        /**
         * This method refines a grid point by generating the children in every dimension
         * of the grid and all their missing ancestors by calling create_gridpoint().
         *
         * @param storage hashmap that stores the gridpoints
         * @param refine_index The index in the hashmap of the point that should be refined
         */
        void refineGridpoint(GridStorage* storage, size_t refine_index);
        /**
         * This method creates a new point on the grid. It checks if some parents or
         * children are needed in other dimensions.
         *
         * @param storage hashmap that stores the gridpoints
         * @param index The point that should be inserted
         */
        void createGridpoint(GridStorage* storage, index_type& index);

        /**
         * Examines the grid points and stores the indices those that can be refined
         * and have maximal indicator values.
         *
         * @param storage hashmap that stores the grid points
         * @param functor a RefinementFunctor specifying the refinement criteria
         * @param refinements_num number of points to refine
         * @param max_indices the array where the point indices should be stored
         * @param max_values the array where the corresponding indicator values
         * should be stored
         */
        virtual void collectRefinablePoints(GridStorage* storage,
                                            RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
                                            RefinementFunctor::value_type* max_values);


        /**
         * Refines the collection of points.
         *
         * @param storage hashmap that stores the grid points
         * @param functor a RefinementFunctor specifying the refinement criteria
         * @param refinements_num number of points to refine
         * @param max_indices the array with the indices of points that should be refined
         * @param max_values the array with the corresponding indicator values
         */
        virtual void refineGridpointsCollection(GridStorage* storage,
                                                RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
                                                RefinementFunctor::value_type* max_values);

    };

  }
}

#endif /* HASHREFINEMENT_HPP */
