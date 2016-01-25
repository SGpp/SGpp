// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef STANDARDGRIDGENERATOR_HPP
#define STANDARDGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * GridGenerator for standard grids without boundaries
     */
    class StandardGridGenerator : public GridGenerator {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's storage object
         */
        StandardGridGenerator(GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~StandardGridGenerator();

        virtual void regular(size_t level);
        virtual void cliques(size_t level, size_t clique_size);
        virtual void full(size_t level);
        virtual void refine(RefinementFunctor* func);
        virtual size_t getNumberOfRefinablePoints();

        virtual void coarsen(CoarseningFunctor* func, DataVector* alpha);
        virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly);
        virtual size_t getNumberOfRemovablePoints();

        virtual void refineMaxLevel(RefinementFunctor* func, size_t maxLevel);
        virtual size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel);

      protected:
        /// pointer to the storage object
        GridStorage* storage;
    };

  }
}

#endif /* STANDARDGRIDGEMERATOR_HPP */
