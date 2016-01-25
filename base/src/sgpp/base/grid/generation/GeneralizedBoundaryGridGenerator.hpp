// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GENERALIZEDTRUNCATEDBOUNDARYGRIDGENERATOR_HPP_
#define GENERALIZEDTRUNCATEDBOUNDARYGRIDGENERATOR_HPP_
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class provides the interface for the grid generation
     * for grids with boundaries, pentagon cut through sub space scheme
     */
    class GeneralizedBoundaryGridGenerator : public GridGenerator {
      public:
        /**
         * Constructor
         *
         * @param storage template type that holds the grid points
         */
        GeneralizedBoundaryGridGenerator(GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~GeneralizedBoundaryGridGenerator();
        /**
         * Creates a regular truncated boundary grid with given level and l_user=1
         * Is the same as the regular truncated grid
         * */
        virtual void regular(size_t level);
        virtual void cliques(size_t level, size_t clique_size);
        virtual void full(size_t level) {};
        /**
         * Creates a super truncated boundary grid with given level and l_user
         * @param level the maximum level of the grid
         * @param l_user the number of fullgrids cut off from the boundaries.
         * */
        virtual void truncated(size_t level, size_t l_user);
        virtual void refine(RefinementFunctor* func) {};
        virtual size_t getNumberOfRefinablePoints() {
          return 0;
        };

        virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) {};
        virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) {};
        virtual size_t getNumberOfRemovablePoints() {
          return 0;
        };

        virtual void refineMaxLevel(RefinementFunctor* func, size_t maxLevel) {
          throw generation_exception("refineMaxLevel is not implemented");
        };
        virtual size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) {
          throw generation_exception("getNumberOfRefinablePointsToMaxLevel is not implemented");
          return 0;
        };

      protected:
        /// Pointer to the grid's storage object
        GridStorage* storage;
    };

  }
}

#endif /* GENERALIZEDTRUNCATEDBOUNDARYGRIDGENERATOR_HPP_ */
