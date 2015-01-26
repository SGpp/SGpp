/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <cmath>

namespace sg {
  namespace base {

    SurplusVolumeRefinementFunctor::SurplusVolumeRefinementFunctor(DataVector* alpha, size_t refinements_num, double threshold) : alpha(alpha), refinements_num(refinements_num), threshold(threshold) {
    }


    SurplusVolumeRefinementFunctor::~SurplusVolumeRefinementFunctor() {
    }

    double SurplusVolumeRefinementFunctor::operator()(GridStorage* storage, size_t seq) {
      return pow(2, static_cast<double>(-((int)storage->get(seq)->getLevelSum()))) * fabs(alpha->get(seq));
    }

    double SurplusVolumeRefinementFunctor::start() {
      return 0.0;
    }

    size_t SurplusVolumeRefinementFunctor::getRefinementsNum() {
      return this->refinements_num;
    }

    double SurplusVolumeRefinementFunctor::getRefinementThreshold() {
      return this->threshold;
    }

  }
}
