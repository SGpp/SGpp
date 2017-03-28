// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULTIPLECLASSREFINEMENTFUNCTOR_HPP
#define MULTIPLECLASSREFINEMENTFUNCTOR_HPP

#include "ZeroCrossingRefinementFunctor.hpp"


#include <vector>
#include <tuple>
#include <algorithm>
#include <sgpp/base/tools/MultipleClassPoint.hpp>

namespace sgpp {
namespace datadriven {

class MultipleClassRefinementFunctor: public ZeroCrossingRefinementFunctor {
 public:
  MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                size_t refinements_num,
                                bool level_penalize,
                                bool pre_compute,
                                double thresh);

  double operator()(base::GridStorage& storage,
                    size_t seq) const override;

  base::Grid* getCombinedGrid();

  void printPointsPlott();
  void printPointsInfo();
  void printScores();
  
  void refine(size_t partCombined);

 private:
  std::vector<sgpp::base::MultipleClassPoint> points;
  base::Grid* multigrid;
  // print different scores
  mutable std::vector<std::string> scoresToPrint;
  bool refineMulti = false;

  void prepareGrid();

  void findCrossings(size_t leftP, size_t rightP, size_t seq, size_t d);
  bool hasChild(const base::HashGridPoint& gp, size_t d, bool left) const;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* MULTIPLECLASSREFINEMENTFUNCTOR_HPP */
