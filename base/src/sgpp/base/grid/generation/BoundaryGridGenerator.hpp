// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TRUNCATEDBOUNDARYGRIDGENERATOR_HPP
#define TRUNCATEDBOUNDARYGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * This class provides the interface for the grid generation
     * for grids with boundaries, pentagon cut through sub space scheme
     */
    class BoundaryGridGenerator : public GridGenerator {
      public:
        /**
         * Constructor
         *
         * @param storage       template type that holds the grid points
         * @param boundaryLevel 1 + how much levels the boundary is coarser than
         *                      the main axes, 0 means one level finer,
         *                      1 means same level,
         *                      2 means one level coarser, etc.
         */
        BoundaryGridGenerator(GridStorage* storage, level_t boundaryLevel = 1);

        /**
         * Destructor
         */
        virtual ~BoundaryGridGenerator() override;

        virtual void regular(size_t level) override;
        virtual void cliques(size_t level, size_t clique_size) override;
        virtual void full(size_t level) override;
        virtual void refine(RefinementFunctor* func) override;
        virtual size_t getNumberOfRefinablePoints() override;

        virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) override;
        virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) override;
        virtual size_t getNumberOfRemovablePoints() override;

        virtual void refineMaxLevel(RefinementFunctor* func, size_t maxLevel) override;
        virtual size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override;

      protected:
        /// Pointer to the grid's storage object
        GridStorage* storage;
        /// 1 + how much levels the boundary is coarser than the main axes
        level_t boundaryLevel;
    };

  }
}

#endif /* TRUNCATEDBOUNDARYGRIDGENERATOR_HPP */
