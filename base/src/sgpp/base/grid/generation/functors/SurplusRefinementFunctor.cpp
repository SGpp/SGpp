/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

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
