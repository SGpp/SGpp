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
 * The zero-crossings are not determinated by pairwise comparing the signes
 * of PDFS f_1 - f_2, but by comparing the dominating class at geometric 
 * neighbors, determined for each dimension separately.
 *
 * Finer levels are pealized by 2^{-|l|}.
 */
class MultipleClassRefinementFunctor: public ZeroCrossingRefinementFunctor {
 public:
  /**
   * Constructor.
   *
   * @param grids Vector of grids. current_grid_index specifies the grid to be refined
   * @param alphas Vector of surpluses related to the grids
   * @param refinements_num Maximum number of refinements done
   * @param partCombined Number of refinement done in the combined grid
   * @param threshold Threshold for refinement scores
   */
  MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                size_t refinements_num,
                                size_t partCombined,
                                double thresh);

  double operator()(base::GridStorage& storage,
                    size_t seq) const override;

  double getTopPercent();
  void setTopPercent(double newPercent);
  double getBorderPenalty();
  void setBorderPenalty(double newPenalty);

  /**
   * Organizes the refinement of the classes
   * uses set parameter to execute the refinement step.
   *
   * Creates the combined grid and starts the refinement.
   */
  void refine();

 private:
  std::vector<sgpp::base::MultipleClassPoint> points;
  base::Grid* multigrid;
  size_t partCombined;
  double topPercent;
  double borderPenalty;

  bool refineMulti = false;
  mutable double borderSum;
  mutable double borderCnt;

  void prepareGrid();

  void findCrossings(size_t leftP, size_t rightP, size_t seq, size_t d);
  bool hasChild(const base::HashGridPoint& gp, size_t d, bool left) const;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* MULTIPLECLASSREFINEMENTFUNCTOR_HPP */
