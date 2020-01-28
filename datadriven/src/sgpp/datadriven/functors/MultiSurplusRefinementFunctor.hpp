// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#ifndef MULTISURPLUSREFINEMENTFUNCTOR_HPP
#define MULTISURPLUSREFINEMENTFUNCTOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>

#include <vector>


namespace sgpp {
namespace datadriven {

/**
 * Wrapper of SurplusRefinementFunctor for multi grid scenarios.
 */
class MultiSurplusRefinementFunctor : public MultiGridRefinementFunctor {
 public:
  /**
   * Constructor.
   *
   * @param grids Vector of grids. current_grid_index specifies the grid to be refined
   * @param alphas Vector of surpluses related to the grids
   * @param refinements_num Maximum number of refinements done
   * @param level_penalize If a level penalizing is multiplied to the score. Here, it determines whether surplus or volume refinement is used
   * @param threshold Threshold for refinement scores
   */
  MultiSurplusRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                size_t refinements_num = 1,
                                bool level_penalize = false,
                                double threshold = 0.0);

  double operator()(base::GridStorage& storage,
                    size_t seq) const override;
  double start() const override;
  size_t getRefinementsNum() const override;
  double getRefinementThreshold() const override;
  ~MultiSurplusRefinementFunctor() override {}

  void setGridIndex(size_t grid_index) override;
  size_t getNumGrids() override;

 protected:
  std::vector<base::Grid*> grids;
  std::vector<base::DataVector*> alphas;

  /**
   * One surplus refinement functor per grid
   */
  std::vector<base::SurplusRefinementFunctor> spFunctors;

  /**
   * One volume refinement functor per grid
   */
  std::vector<base::SurplusVolumeRefinementFunctor> spvFunctors;

  size_t current_grid_index;
  bool level_penalize;
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* ZEROCROSSINGFUNCTOR_HPP */
