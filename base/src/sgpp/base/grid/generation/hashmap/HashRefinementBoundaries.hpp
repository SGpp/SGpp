/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHREFINEMENTBOUNDARIES_HPP
#define HASHREFINEMENTBOUNDARIES_HPP

#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Standard free refinement class for sparse grids with boundaries
     */
    class HashRefinementBoundaries: public AbstractRefinement {
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

        void refineGridpoint1D(GridStorage* storage, AbstractRefinement::index_type& index, size_t d);

      protected:
        /**
         * This method refines a grid point be generating the children in every dimension
         * of the grid.
         *
         * @param storage hashmap that stores the gridpoints
         * @param refine_index the index in the hashmap of the point that should be refined
         */
        void refineGridpoint(GridStorage* storage, size_t refine_index);


        /**
         * Wrapper for the two functions create_gridpoint_general and
         * create_gridpoint_levelZeroConsistency which have both to be
         * executed if a gridpoint is refined
         *
         * @param storage hashmap that stores the gridpoinrs
         * @param index the point that should be inserted
         */
        void createGridpoint(GridStorage* storage, index_type& index);


        /**
         * This method creates a new point on the grid. It checks if some parents or
         * children are needed in other dimensions.
         *
         * @param storage hashmap that stores the gridpoinrs
         * @param index the point that should be inserted
         */
        void createGridpointGeneral(GridStorage* storage, index_type& index);


        /**
         * Assures that we have always a consistent grid with both functions
         * 0,0 and 0,1 on level zero
         *
         * @param storage hashmap that stores the gridpoinrs
         * @param index the point that should be inserted
         */
        void createGridpointLevelZeroConsistency(GridStorage* storage, index_type& index);

        /**
               * Creates children grid points along single direction
               *
               * @param index The point that should be refined
               * @param d direction
               * @param storage hashmap that stores the gridpoints
               * @param source_index index value in the dimension d
               * @param source_level level value in the dimension d
               */
        void createGridpoint1D(index_type& index,
                               size_t d, GridStorage* storage,
                               index_t& source_index, level_t& source_level);

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

#endif /* HASHREFINEMENTBOUNDARIES_HPP */
