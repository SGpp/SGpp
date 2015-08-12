// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REFINEMENTFUNCTOR_HPP
#define REFINEMENTFUNCTOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Abstract class that defines the interface that refinement functors have to provide
     */
    class RefinementFunctor {
      public:
        typedef float_t value_type;

        /**
         * Constructor
         */
        RefinementFunctor() {}

        /**
         * Destructor
         */
        virtual ~RefinementFunctor() {}

        /**
         * This should be returning a refinement value for every grid point.
         * The point with the highest value will be refined first.
         *
         * @param storage pointer to the grids storage object
         * @param seq sequence number in the coefficients array
         *
         * @return refinement value
         */
        virtual float_t operator()(GridStorage* storage, size_t seq) = 0;

        /**
         * Returns the lower bound of refinement criterion (e.g., alpha or error) (lower bound).
         * The refinement value of grid points to be refined have to be larger than this value
         *
         * @return lower bound
         */
        virtual float_t start() = 0;

        /**
         * Returns the maximal number of points that should be refined.
         *
         * The maximal number of points to refine is set in the constructor of the implementing class.
         *
         * @return number of points that should refined. Default value: 1.
         */
        virtual size_t getRefinementsNum() {
          return 1;
        }

        /**
         * Returns the threshold for refinement.
         *
         * Only the grid points with absolute value of refinement criterion (e.g., alpha or error) greater
         * or equal to this threshold will be refined.
         *
         * @return threshold value for refinement. Default value: 0.
         */
        virtual float_t getRefinementThreshold() = 0;

        /**
         * Returns the total sum of local (error) indicators used for refinement
         *
         * @param storage pointer to the grids storage object
         * @return total sum of local (error) indicators used for refinement
         */
        virtual float_t getTotalRefinementValue(GridStorage* storage) {
          float_t sum = 0;
          GridStorage::grid_map_iterator end_iter = storage->end();

          for (GridStorage::grid_map_iterator iter = storage->begin();
               iter != end_iter; iter++) {
            sum += operator()(storage, iter->second);
          }

          return sum;
        }
    };

  }
}

#endif /* REFINEMENTFUNCTOR_HPP */