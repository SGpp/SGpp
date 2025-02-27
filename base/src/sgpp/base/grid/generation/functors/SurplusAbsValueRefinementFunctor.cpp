// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusAbsValueRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace base {

SurplusAbsValueRefinementFunctor::SurplusAbsValueRefinementFunctor(Grid& grid, DataVector& alpha,
                                                                   size_t refinements_num,
                                                                   double threshold)
    : grid(grid),
      alpha(alpha),
      evals(grid.getSize()),
      refinements_num(refinements_num),
      threshold(threshold) {
  // collect coordinates of all grid points
  base::DataMatrix gridPoints(grid.getSize(), grid.getDimension());
  base::DataVector p(grid.getDimension());

  auto& gridStorage = grid.getStorage();
  for (size_t i = 0; i < grid.getSize(); ++i) {
    gridStorage.getPoint(i).getStandardCoordinates(p);
    gridPoints.setRow(i, p);
  }

  op_factory::createOperationMultipleEvalDefault(grid, gridPoints)->eval(alpha, evals);
}

SurplusAbsValueRefinementFunctor::~SurplusAbsValueRefinementFunctor() {}

double SurplusAbsValueRefinementFunctor::operator()(GridStorage& storage, size_t seq) const {
  return fabs(evals[seq]) * fabs(alpha[seq]);
}

double SurplusAbsValueRefinementFunctor::start() const { return 0.0; }

size_t SurplusAbsValueRefinementFunctor::getRefinementsNum() const { return this->refinements_num; }

double SurplusAbsValueRefinementFunctor::getRefinementThreshold() const { return this->threshold; }

}  // namespace base
}  // namespace sgpp
