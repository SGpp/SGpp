// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#ifndef MULTIGRIDCOARSENINGFUNCTOR_HPP
#define MULTIGRIDCOARSENINGFUNCTOR_HPP

#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>


namespace sgpp {
    namespace datadriven {

/**
 * Abstract super-class for coarsening functors operating on multiple
 * grids.
 */
        class MultiGridCoarseningFunctor : public base::CoarseningFunctor {
        public:
            /**
             * Sets the index (into the vector of grids) of the grid to be coarsened
             *
             * @param grid_index The index of the grid to be coarsened
             */
            virtual void setGridIndex(size_t grid_index) = 0;


            /**
             * Returns the number of grids the functor can / does coarsen
             */
            virtual size_t getNumGrids() = 0;

            /**
             * Used if expensive computations (eg. grid evaluations)
             * are cached, usually for one coarsen step.
             */
            virtual void preComputeEvaluations() { }

            virtual ~MultiGridCoarseningFunctor() { }
        };
    }  // namespace datadriven
}  // namespace sgpp

#endif

