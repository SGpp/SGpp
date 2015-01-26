/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Aliz Nagy

#ifndef TRUNCATEDTRAPEZOIDGRIDGENERATOR_HPP_
#define TRUNCATEDTRAPEZOIDGRIDGENERATOR_HPP_
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

namespace sg {
  namespace base {

    /**
     * This class provides the interface for the grid generation
     * for grids with boundaries, pentagon cut through sub space scheme
     */
    class TruncatedTrapezoidGridGenerator : public GridGenerator {
      public:
        /**
         * Constructor
         *
         * @param storage template type that holds the grid points
         */
        TruncatedTrapezoidGridGenerator(GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~TruncatedTrapezoidGridGenerator();
        /**
         * Creates a regular trapezoid boundary grid with given level and l_user=1
         * Is the same as the regular trapezoid grid
         * */
        virtual void regular(int level);
        virtual void cliques(int level, size_t clique_size);
        virtual void full(int level) {};
        /**
         * Creates a super trapezoid boundary grid with given level and l_user
         * @param level the maximum level of the grid
         * @param l_user the number of fullgrids cut off from the boundaries.
         * */
        virtual void truncated(int level, int l_user);
        virtual void refine(RefinementFunctor* func) {};
        virtual size_t getNumberOfRefinablePoints() {
          return 0;
        };

        virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) {};
        virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) {};
        virtual size_t getNumberOfRemovablePoints() {
          return 0;
        };

        virtual void refineMaxLevel(RefinementFunctor* func, int maxLevel) {
          throw generation_exception("refineMaxLevel is not implemented");
        };
        virtual size_t getNumberOfRefinablePointsToMaxLevel(int maxLevel) {
          throw generation_exception("getNumberOfRefinablePointsToMaxLevel is not implemented");
          return 0;
        };

      protected:
        /// Pointer to the grid's storage object
        GridStorage* storage;
    };

  }
}

#endif /* TRUNCATEDTRAPEZOIDGRIDGENERATOR_HPP_ */
