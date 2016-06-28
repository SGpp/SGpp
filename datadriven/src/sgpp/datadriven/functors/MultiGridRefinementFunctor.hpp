// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#ifndef MULTIGRIDREFINMENTFUNCTOR_HPP
#define MULTIGRIDREFINMENTFUNCTOR_HPP

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>


namespace sgpp {
namespace datadriven {
  /**
   * Abstract super-class for refinement functors operating on multiple
   * grids.
   */
  class MultiGridRefinementFunctor : public base::RefinementFunctor {
    public:
      /**
       * Sets the index (into the vector of grids) of the grid to be refined
       *
       * @param grid_index The index of the grid to be refined
       */
      virtual void setGridIndex(size_t grid_index) = 0;

      virtual size_t getNumGrids() = 0;

      /**
       * Used if expensive computations (eg. grid evaluations)
       * are cached, usually for one refinement step.
       */
      virtual void preComputeEvaluations() { }

      virtual ~MultiGridRefinementFunctor() { }
  };
} // namespace datadriven
} // namespace sgpp


#endif
