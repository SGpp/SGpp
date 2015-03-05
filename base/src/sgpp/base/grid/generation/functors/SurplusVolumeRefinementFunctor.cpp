// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <cmath>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    SurplusVolumeRefinementFunctor::SurplusVolumeRefinementFunctor(DataVector* alpha, size_t refinements_num, float_t threshold) : alpha(alpha), refinements_num(refinements_num), threshold(threshold) {
    }


    SurplusVolumeRefinementFunctor::~SurplusVolumeRefinementFunctor() {
    }

    float_t SurplusVolumeRefinementFunctor::operator()(GridStorage* storage, size_t seq) {
      return pow(2, static_cast<float_t>(-((int)storage->get(seq)->getLevelSum()))) * fabs(alpha->get(seq));
    }

    float_t SurplusVolumeRefinementFunctor::start() {
      return 0.0;
    }

    size_t SurplusVolumeRefinementFunctor::getRefinementsNum() {
      return this->refinements_num;
    }

    float_t SurplusVolumeRefinementFunctor::getRefinementThreshold() {
      return this->threshold;
    }

  }
}