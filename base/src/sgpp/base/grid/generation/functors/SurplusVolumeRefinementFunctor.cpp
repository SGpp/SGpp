// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>


namespace sgpp {
namespace base {

SurplusVolumeRefinementFunctor::SurplusVolumeRefinementFunctor(
  DataVector& alpha, size_t refinements_num, double threshold) : alpha(alpha),
  refinements_num(refinements_num), threshold(threshold) {
}


SurplusVolumeRefinementFunctor::~SurplusVolumeRefinementFunctor() {
}

double SurplusVolumeRefinementFunctor::operator()(GridStorage& storage,
    size_t seq) const {
  return pow(2, static_cast<double>(
               -(static_cast<int>(storage.get(seq)->getLevelSum())))) * fabs(alpha[seq]);
}

double SurplusVolumeRefinementFunctor::start() const {
  return 0.0;
}

size_t SurplusVolumeRefinementFunctor::getRefinementsNum() const {
  return this->refinements_num;
}

double SurplusVolumeRefinementFunctor::getRefinementThreshold() const {
  return this->threshold;
}

}  // namespace base
}  // namespace sgpp
