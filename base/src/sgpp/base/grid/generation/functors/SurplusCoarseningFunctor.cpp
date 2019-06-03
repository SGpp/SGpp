// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

SurplusCoarseningFunctor::SurplusCoarseningFunctor(DataVector& alpha,
    size_t removements_num, double threshold) : alpha(alpha),
  removements_num(removements_num), threshold(threshold) {
}

SurplusCoarseningFunctor::~SurplusCoarseningFunctor() {}


double SurplusCoarseningFunctor::operator()(GridStorage& storage, size_t seq) const {
  return fabs(alpha[seq]);
}

double SurplusCoarseningFunctor::start() const {
  return 1.0;
}

size_t SurplusCoarseningFunctor::getRemovementsNum() const {
  return this->removements_num;
}

double SurplusCoarseningFunctor::getCoarseningThreshold() const {
  return this->threshold;
}

}  // namespace base
}  // namespace sgpp
