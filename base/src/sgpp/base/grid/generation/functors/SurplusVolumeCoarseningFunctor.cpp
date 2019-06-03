// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusVolumeCoarseningFunctor.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>


namespace sgpp {
namespace base {

SurplusVolumeCoarseningFunctor::SurplusVolumeCoarseningFunctor(DataVector& alpha,
    size_t removements_num, double threshold) : alpha(alpha),
  removements_num(removements_num), threshold(threshold) {
}

SurplusVolumeCoarseningFunctor::~SurplusVolumeCoarseningFunctor() {}


double SurplusVolumeCoarseningFunctor::operator()(GridStorage& storage, size_t seq) const {
  return pow(2, static_cast<double>(
               -(static_cast<int>(storage.getPoint(seq).getLevelSum())))) * fabs(alpha[seq]);
}

double SurplusVolumeCoarseningFunctor::start() const {
  return 1.0;
}

size_t SurplusVolumeCoarseningFunctor::getRemovementsNum() const {
  return this->removements_num;
}

double SurplusVolumeCoarseningFunctor::getCoarseningThreshold() const {
  return this->threshold;
}

}  // namespace base
}  // namespace sgpp
