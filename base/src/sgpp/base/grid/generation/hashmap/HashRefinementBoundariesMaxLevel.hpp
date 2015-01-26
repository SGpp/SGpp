/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASHREFINEMENTBOUNDARIESMAXLEVEL_HPP
#define HASHREFINEMENTBOUNDARIESMAXLEVEL_HPP

#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

namespace sg {
  namespace base {

    /**
     * Standard free refinement class for sparse grids with boundaries
     * with a maximal level depth of refinement
     */
    class HashRefinementBoundariesMaxLevel: public HashRefinementBoundaries {
      public:

        /**
         * Performs the refinement on grid
         *
         * @param storage hashmap that stores the grid points
         * @param functor a function used to determine if refinement is needed
         * @param maxLevel no points on higher levels than maxLevel will be created
         */
        void refineToMaxLevel(GridStorage* storage, RefinementFunctor* functor, unsigned int maxLevel);


        /**
         * Calculates the number of points, which can be refined
         *
         * @param storage hashmap that stores the grid points
         * @param maxLevel no points on higher levels than maxLevel will be created
         */
        size_t getNumberOfRefinablePointsToMaxLevel(GridStorage* storage, unsigned int maxLevel);


        /**
         * Refine one grid point along a single direction
         * @param storage hashmap that stores the grid points
         * @param index point to refine
         * @param d direction
         * @param maxLevel no points on higher levels than maxLevel will be created
               */
        void refineGridpoint1D(GridStorage* storage, index_type& index,
                               size_t d, unsigned int maxLevel);


      protected:
        /**
         * This method refines a grid point be generating the children in every dimension
         * of the grid.
         *
         * @param storage hashmap that stores the gridpoints
         * @param refine_index the index in the hashmap of the point that should be refined
         * @param maxLevel no points on higher levels than maxLevel will be created
         */
        void refineGridpoint(GridStorage* storage, size_t refine_index, unsigned int maxLevel);


    };

  }
}

#endif /* HASHREFINEMENTBOUNDARIESMAXLEVEL_HPP */
