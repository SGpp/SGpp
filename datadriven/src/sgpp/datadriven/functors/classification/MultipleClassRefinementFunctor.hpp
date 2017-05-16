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

class MultipleClassRefinementFunctor: public ZeroCrossingRefinementFunctor {
 public:
  MultipleClassRefinementFunctor(std::vector<base::Grid*> grids,
                                std::vector<base::DataVector*> alphas,
                                size_t refinements_num,
                                double thresh);

  double operator()(base::GridStorage& storage,
                    size_t seq) const override;

  double getTopPercent();
  void setTopPercent(double newPercent);
  double getBorderPenalty();
  void setBorderPenalty(double newPenalty);

  void printPointsPlott();
  void printPointsInfo();

  void refine(size_t partCombined);

 private:
  std::vector<sgpp::base::MultipleClassPoint> points;
  base::Grid* multigrid;
  bool refineMulti = false;
  mutable double borderSum;
  mutable double borderCnt;

  void prepareGrid();

  void findCrossings(size_t leftP, size_t rightP, size_t seq, size_t d);
  bool hasChild(const base::HashGridPoint& gp, size_t d, bool left) const;

  double topPercent;
  double borderPenalty;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* MULTIPLECLASSREFINEMENTFUNCTOR_HPP */
