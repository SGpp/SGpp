// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef SQUAREROOTGRIDGENERATOR_HPP_
#define SQUAREROOTGRIDGENERATOR_HPP_
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class provides the interface for the grid generation
     * for grids with boundaries, pentagon cut through sub space scheme
     */
    class SquareRootGridGenerator : public GridGenerator {
      public:
        /**
         * Constructor
         *
         * @param storage template type that holds the grid points
         */
        SquareRootGridGenerator(GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~SquareRootGridGenerator();

        virtual void regular(int level);
        virtual void cliques(int level, size_t clique_size);
        virtual void full(int level) {};
        virtual void refine(RefinementFunctor* func) {};
        virtual size_t getNumberOfRefinablePoints() {
          return 0;
        };

        virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) {};
        virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) {};
        virtual size_t getNumberOfRemovablePoints() {
          return 0;
        };

        virtual void refineMaxLevel(RefinementFunctor* func, int maxLevel) {};
        virtual size_t getNumberOfRefinablePointsToMaxLevel(int maxLevel) {
          return 0;
        };

      protected:
        /// Pointer to the grid's storage object
        GridStorage* storage;
    };

  }
}

#endif /* SQUAREROOTGRIDGENERATOR_HPP_ */