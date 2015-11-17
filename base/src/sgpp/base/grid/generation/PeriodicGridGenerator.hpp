// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PERIODICGRIDGENERATOR_HPP
#define PERIODICGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * GridGenerator for periodic grids with boundaries
     */
    class PeriodicGridGenerator : public GridGenerator {
      public:
        /**
         * Constructor
         *
         * @param storage the grid's storage object
         */
        PeriodicGridGenerator(GridStorage* storage);

        /**
         * Destructor
         */
        virtual ~PeriodicGridGenerator() override;

        virtual void regular(size_t level) override;
        virtual void full(size_t level) override;
        virtual void refine(RefinementFunctor* func) override;
        virtual void cliques(size_t level, size_t clique_size) override;
        virtual size_t getNumberOfRefinablePoints() override;

        virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) override;
        virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) override;
        virtual size_t getNumberOfRemovablePoints() override;

        virtual void refineMaxLevel(RefinementFunctor* func, size_t maxLevel) override;
        virtual size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override;

      protected:
        /// pointer to the storage object
        GridStorage* storage;
    };

  }
}


#endif /* PERIODICGRIDGENERATOR_HPP */
