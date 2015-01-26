// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    SurplusRefinementFunctor::SurplusRefinementFunctor(DataVector* alpha, size_t refinements_num, double threshold) : alpha(alpha), refinements_num(refinements_num), threshold(threshold) {
    }


    SurplusRefinementFunctor::~SurplusRefinementFunctor() {
    }

    double SurplusRefinementFunctor::operator()(GridStorage* storage, size_t seq) {
      double val = fabs(alpha->get(seq));
      // std::cout << seq << ", ";
      return val;
    }

    double SurplusRefinementFunctor::start() {
      return 0.0;
    }

    size_t SurplusRefinementFunctor::getRefinementsNum() {
      return this->refinements_num;
    }

    double SurplusRefinementFunctor::getRefinementThreshold() {
      return this->threshold;
    }

  }
}