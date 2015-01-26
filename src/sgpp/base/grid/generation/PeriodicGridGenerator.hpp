/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#ifndef PERIODICGRIDGENERATOR_HPP
#define PERIODICGRIDGENERATOR_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/GridGenerator.hpp"

namespace sg {
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
        virtual ~PeriodicGridGenerator();

        virtual void regular(int level);
        virtual void full(int level);
        virtual void refine(RefinementFunctor* func);
        virtual void cliques(int level, size_t clique_size);
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


#endif /* PERIODICGRIDGENERATOR_HPP */
