// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace sgpp {
namespace datadriven {

  MultiSurplusRefinementFunctor::
  MultiSurplusRefinementFunctor(std::vector<base::Grid*> g,
                                std::vector<base::DataVector*> a,
                                size_t r_num,
                                bool level_penalize,
                                double thresh) :
    grids(g), alphas(a), current_grid_index(-1),
    level_penalize(level_penalize) {
    for (size_t i = 0; i < grids.size(); i++) {
      spFunctors.push_back(base::SurplusRefinementFunctor(*alphas.at(i),
                                                          r_num,
                                                          thresh));
      spvFunctors.push_back(base::
                            SurplusVolumeRefinementFunctor(*alphas.at(i),
                                                           r_num,
                                                           thresh));
    }
  }

  double MultiSurplusRefinementFunctor::operator()(base::GridStorage&
                                                   storage,
                                                   size_t seq) const {
    size_t cgi = current_grid_index;
    if (level_penalize) {
      return spvFunctors.at(cgi)(storage, seq);
    } else {
      return spFunctors.at(cgi)(storage, seq);
    }
  }

  double MultiSurplusRefinementFunctor::start() const {
    if (level_penalize) {
      return spvFunctors.at(current_grid_index).start();
    } else {
      return spFunctors.at(current_grid_index).start();
    }
  }

  size_t MultiSurplusRefinementFunctor::getRefinementsNum() const {
    if (level_penalize) {
      return spvFunctors.at(current_grid_index).getRefinementsNum();
    } else {
      return spFunctors.at(current_grid_index).getRefinementsNum();
    }
  }

  double MultiSurplusRefinementFunctor::getRefinementThreshold() const {
    if (level_penalize) {
      return spvFunctors.at(current_grid_index).getRefinementThreshold();
    } else {
      return spFunctors.at(current_grid_index).getRefinementThreshold();
    }
  }

  void MultiSurplusRefinementFunctor::setGridIndex(size_t grid_index) {
    this->current_grid_index = grid_index;
  }

  size_t MultiSurplusRefinementFunctor::getNumGrids() {
    return this->grids.size();
  }

} // namespace datadriven
} // namespace sgpp
