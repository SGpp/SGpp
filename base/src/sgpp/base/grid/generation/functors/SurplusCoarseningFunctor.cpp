// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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