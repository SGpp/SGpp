// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <map>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Classification refinement builds up the combined grid of all classes,
 * then iterates through it always keeping track of the neighbors to
 * detect class changes between neighbors. Only the children in the
 * relevant dimension and direction will be added to the grid.
 */
class ClassificationRefinementFunctor : public MultiGridRefinementFunctor {
 public:
  /**
   * Constructor.
   *
   * @param grids Vector of grids. current_grid_index specifies the grid to be refined
   * @param alphas Vector of surpluses related to the grids
   * @param priors Vector of priors related to the classificator
   * @param refinements_num Maximum number of refinements done
   * @param coarsenings_num Maximum number of coarsenings done
   * @param level_penalize If a level penalizing is multiplied to the score (2^{|l|_1})
   * @param thresholdType type of the threshold: absolute score value or relative score value
   * @param refinementThreshold Threshold for refinement scores
   * @param coarseningThreshold Threshold for coarsening scores
   * @param coarsenInitialPoints Whether or not to coarsen initial grid points
   * @param minimumCoarseningIndex Mimimum index of grid points that can be coarsened, used when
   * coarsenInitialPoints = false.
   */
  ClassificationRefinementFunctor(std::vector<base::Grid*> grids,
                                  std::vector<base::DataVector*> alphas, std::vector<double> priors,
                                  size_t refinements_num = 1, size_t coarsenings_num = 1,
                                  bool level_penalize = true,
                                  sgpp::base::AdaptivityThresholdType thresholdType =
                                      sgpp::base::AdaptivityThresholdType::Absolute,
                                  double refinementThreshold = 0.0,
                                  double coarseningThreshold = 1.0,
                                  bool coarsenInitialPoints = true,
                                  size_t minimumCoarseningIndex = 0);

  double operator()(base::GridStorage& storage, size_t seq) const override;
  double start() const override;
  size_t getRefinementsNum() const override;
  double getRefinementThreshold() const override;
  ~ClassificationRefinementFunctor() override {}

  void setGridIndex(size_t grid_index) override;
  size_t getNumGrids() override;
  void preComputeEvaluations() override;

  /**
   * Refine and Coarsen all grids of the model.
   * @returns vector that for each class contains a list of indices of deleted grid points
   */
  std::vector<std::vector<size_t>> adaptAllGrids();

 protected:
  std::vector<base::Grid*> grids;
  std::vector<base::DataVector*> alphas;
  std::vector<double> priors;

  size_t refinements_num;

  /**
   * Maximum number of coarsened points per iteration.
   */
  size_t coarsenings_num;

  bool level_penalize;
  sgpp::base::AdaptivityThresholdType thresholdType;
  double refinementThreshold;
  double coarseningThreshold;

  /**
   * Whether or not to coarsen initial grid points.
   */
  bool coarsenInitialPoints;

  /**
   * Maximum index of the initial grid points.
   */
  size_t minimumCoarseningIndex;

  base::GridStorage total_grid;

  std::vector<std::tuple<size_t, size_t, size_t, bool, base::GridPoint::level_type>> neighborRels;

  /**
   * Stores grid evaluations at all grids (vector) at the union
   * of grid point coordinates over all grids (hashed by string representation
   * in the map)
   */
  std::vector<std::map<std::string, double>> pre_comp_evals;

  /**
   * Gets the evaluations of all grids at the coords of seq
   */
  std::vector<double> getEvalVector(size_t ind, size_t seq) const;

  int sgn(double d) const;

  // Utility for finding geometric neighbors

  /**
   * Iterate through the sparse grid tree
   */
  void stepDown(size_t d, size_t minDim, base::HashGridPoint& gp,
                std::vector<std::pair<base::HashGridPoint, base::HashGridPoint>>& neighbors);

  /**
   * Stores the relation of two found neighbors
   */
  void collectNeighbors(base::HashGridPoint leaf, base::HashGridPoint neighbor, size_t dim,
                        bool isLeft);

  /**
   * Checks if a refinement candidate should be inserted into the map of candidates.
   */
  void insertCoarseningCandidate(
      size_t y,
      std::vector<std::multimap<double, std::tuple<size_t, size_t, bool>>>& classMapsCoarsening,
      size_t leafSeqNumber, double score, std::tuple<size_t, size_t, bool> candidate);
};
}  // namespace datadriven
}  // namespace sgpp
