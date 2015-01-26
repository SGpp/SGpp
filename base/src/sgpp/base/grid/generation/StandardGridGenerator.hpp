/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef STANDARDGRIDGENERATOR_HPP
#define STANDARDGRIDGENERATOR_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/GridGenerator.hpp"

namespace sg {
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

        virtual void regular(int level);
        virtual void cliques(int level, size_t clique_size);
        virtual void full(int level);
        virtual void refine(RefinementFunctor* func);
        virtual size_t getNumberOfRefinablePoints();

        virtual void coarsen(CoarseningFunctor* func, DataVector* alpha);
        virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly);
        virtual size_t getNumberOfRemovablePoints();

        virtual void refineMaxLevel(RefinementFunctor* func, int maxLevel);
        virtual size_t getNumberOfRefinablePointsToMaxLevel(int maxLevel);

      protected:
        /// pointer to the storage object
        GridStorage* storage;
    };

  }
}

#endif /* STANDARDGRIDGEMERATOR_HPP */
