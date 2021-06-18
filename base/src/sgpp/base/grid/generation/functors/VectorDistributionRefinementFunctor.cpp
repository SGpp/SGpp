// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/VectorDistributionRefinementFunctor.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

VectorDistributionRefinementFunctor::VectorDistributionRefinementFunctor(
    DataMatrix& alphas, sgpp::base::DistributionsVector pdfs, size_t refinements_num,
    double threshold)
    : alphas(alphas), pdfs(pdfs), refinements_num(refinements_num), threshold(threshold) {
  sgpp::base::DataMatrix bounds = pdfs.getBounds();
  sampleVector = pdfs.sample();
  // map sampleVector to uni intervals to make it comparable with the grid points
  for (size_t d = 0; d < pdfs.getSize(); d++) {
    sampleVector[d] = (sampleVector[d] - bounds.get(d, 0)) / (bounds.get(d, 1) - bounds.get(d, 0));
  }
}

VectorDistributionRefinementFunctor::~VectorDistributionRefinementFunctor() {}

double VectorDistributionRefinementFunctor::operator()(GridStorage& storage, size_t seq) const {
  sgpp::base::DataVector point = storage.getPointCoordinates(seq);
  point.sub(sampleVector);
  double distance = point.l2Norm();
  // the largest value is picked for refinement, so we return the inverse of the distance
  return 1.0 / distance;
}

double VectorDistributionRefinementFunctor::start() const { return 0.0; }

size_t VectorDistributionRefinementFunctor::getRefinementsNum() const {
  return this->refinements_num;
}

double VectorDistributionRefinementFunctor::getRefinementThreshold() const {
  return this->threshold;
}

}  // namespace base
}  // namespace sgpp
