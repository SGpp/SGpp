// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULTIPLECLASSREFINEMENTFUNCTOR_HPP
#define MULTIPLECLASSREFINEMENTFUNCTOR_HPP

#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/base/tools/MultipleClassPoint.hpp>

#include <vector>
#include <tuple>
#include <string>
#include <algorithm>

namespace sgpp {
namespace datadriven {

/**
 * Multiple class refinement is based on the zero-crossing based refinement.
 * The zero-crossings are not determinate by pairwise comparing the sings
 * of PDFS f_1 - f_2, but by comparing the dominating class at geometric 
 * neighbors, determined for each dimension separately.
 *
 * Finer levels are penalized by 2^{-|l|}.
 */
class MultipleClassRefinementFunctor: public ZeroCrossingRefinementFunctor {
 public:
  /**
   * Constructor.
   *
   * @param grids Vector of grids. current_grid_index specifies the grid to be refined
   * @param alphas Vector of surpluses related to the grids
   * @param priors Vector of priors related to the classificator
   * @param refinements_num Maximum number of refinements done
   * @param partCombined Number of refinement done in the combined grid
   * @param thresh Threshold for refinement scores
   */
  MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                std::vector<double> priors,
                                size_t refinements_num,
                                size_t partCombined,
                                double thresh);

  double operator()(base::GridStorage& storage,
                    size_t seq) const override;

  /**
   * Gets the range in which densities are considered to be close.
   * Gives a percentage [0,1].
   *
   * @return Range in which densities are considered close
   */
  double getTopPercent();
  /**
   * Sets the range in which densities are considered to be close.
   * Set the percentage [0,1] the densities need to have, compared to the
   * density of the dominating class.
   *
   * @param newPercent Set the new range
   */
  void setTopPercent(double newPercent);
  /**
   * Gets the factor, the borders are penalized with.
   *
   * @return The border penalty
   */
  double getBorderPenalty();
  /**
   * Sets the factor, used to penalize the borders.
   * The border score is multiplied with the factor before beeing added
   * to the overall score.
   *
   * @param newPenalty Factor to penalize the borders
   */
  void setBorderPenalty(double newPenalty);

  /**
   * Organizes the refinement of the classes
   * uses set parameter to execute the refinement step.
   *
   * Creates the combined grid and starts the refinement.
   */
  void refine();

 private:
  // Vector saving the additional informations for each GridPoint
  std::vector<sgpp::base::MultipleClassPoint> points;
  // Combined grid, containing all points of the grids from each class
  base::Grid* multigrid;
  // Points refined in the combined grid, not in the classes
  size_t partCombined;
  // Range for close densities
  double topPercent;
  // factor to penalize the boundaries
  double borderPenalty;

  // what grid is refined, if false, use current_grid_index
  bool refineMulti = false;
  mutable double borderSum;
  mutable double borderCnt;

  /**
   * Creates the combined grid for the following refinement step.
   * Sets the variable multigrid.
   */
  void prepareGrid();

  /**
   * Finds all neighbors with different dominating class.
   * Sets neighbors/borders in multigrid.
   * Recursively searches all GridPoints and all dimensions
   *
   * @param leftP Sequence number of left neighbor, greater than GridSize if not found
   * @param rightP Sequence number of right neighbor, greater than GridSize if not found
   * @param seq Sequence number of the currently calculated GridPoint
   * @param d Dimension currently searched
   */
  void findCrossings(size_t leftP, size_t rightP, size_t seq, size_t d);
  bool hasChild(const base::HashGridPoint& gp, size_t d, bool left) const;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* MULTIPLECLASSREFINEMENTFUNCTOR_HPP */
