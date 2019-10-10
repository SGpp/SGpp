// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULTIPLECLASSREFINEMENT_HPP
#define MULTIPLECLASSREFINEMENT_HPP

#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/MultipleClassPoint.hpp>

#include <iostream>
#include <tuple>
#include <cmath>
#include <vector>
#include <algorithm>


namespace sgpp {
namespace base {

/**
 * Refinement class for sparse grids.
 * Used by the MultipleClassRefinementFuntor.
 * Scores all grid points and refines points based on the information
 * given in the vector of MulitpleClassPoint
 */
class MultipleClassRefinement : public HashRefinement {
 public:
  /**
   * Constructor.
   *
   * @param grid Combined grid. current_grid_index specifies the grid to be refined
   * @param pts Vector of MultipleClassPoint with additional information
   * @param classGrids Vector of grids
   * @param borderSum Sum of the border scores for all points
   * @param borderCnt amount of points scored towards the border
   * @param topPercent range when densities are concidered close
   */
    MultipleClassRefinement(Grid& grid,
        std::vector<sgpp::base::MultipleClassPoint>* pts,
        std::vector<Grid*>& classGrids,
        double &borderSum, double &borderCnt, double topPercent);
    ~MultipleClassRefinement() override {}

 protected:
  void refineGridpoint(GridStorage& storage, size_t refine_index) override;
  void collectRefinablePoints(GridStorage& storage,
        RefinementFunctor& functor,
        AbstractRefinement::refinement_container_type& collection) override;

 private:
    // Additional data for combined grid
    std::vector<sgpp::base::MultipleClassPoint>* points;
    // Combined grid
    Grid& multigrid;
    // Grids of the classes
    std::vector<Grid*>& grids;
    double &borderSum;
    double &borderCnt;
    // Range for close densities
    double topPercent;

    /**
     * Adds the given Gridpoint into the given GridStorage if not already inside.
     *
     * @param storage GridStorage the Gridpoint should be inserted
     * @param point Gridpoint to insert into given GridStorage
     */
    void addGridpoint(GridStorage& storage, GridPoint& point);
};

} /* namespace base */
} /* namespace sgpp */

#endif /* MULTIPLECLASSREFINEMENT_HPP */
