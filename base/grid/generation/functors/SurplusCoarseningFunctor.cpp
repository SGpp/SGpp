/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/grid/generation/functors/SurplusCoarseningFunctor.hpp"

namespace sg {
  namespace base {

    SurplusCoarseningFunctor::SurplusCoarseningFunctor(DataVector* alpha, size_t removements_num, double threshold) : alpha(alpha), removements_num(removements_num), threshold(threshold) {
    }

    SurplusCoarseningFunctor::~SurplusCoarseningFunctor() {}


    double SurplusCoarseningFunctor::operator()(GridStorage* storage, size_t seq) {
      return fabs(alpha->get(seq));
    }

    double SurplusCoarseningFunctor::start() {
      return 1.0;
    }

    size_t SurplusCoarseningFunctor::getRemovementsNum() {
      return this->removements_num;
    }

    double SurplusCoarseningFunctor::getCoarseningThreshold() {
      return this->threshold;
    }

  }
}
