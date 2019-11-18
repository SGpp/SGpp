// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>

#include <string>
#include <vector>
#include <list>

using sgpp::base::RefinementFunctor;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::SurplusVolumeRefinementFunctor;
using sgpp::base::RefinementFunctorType;

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimation::ModelFittingDensityEstimation() : refinementsPerformed{0} {}

RefinementFunctor *ModelFittingDensityEstimation::getRefinementFunctor() {
  sgpp::base::AdaptivityConfiguration& refinementConfig = this->config->getRefinementConfig();
  switch (refinementConfig.refinementFunctorType) {
    case RefinementFunctorType::Surplus : {
      return new SurplusRefinementFunctor(alpha, config->getRefinementConfig().noPoints_,
                                             config->getRefinementConfig().threshold_);
    }
    case RefinementFunctorType::SurplusVolume : {
      return new SurplusVolumeRefinementFunctor(alpha, config->getRefinementConfig().noPoints_,
          config->getRefinementConfig().threshold_);
    }
    case RefinementFunctorType::DataBased : {
      std::string errorMessage = "Unsupported refinement functor type DataBased "
          "for density estimation!";
      throw new application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::ZeroCrossing : {
      std::string errorMessage = "Unsupported refinement functor type ZeroCrossing "
          "for density estimation!";
      throw new application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::MultipleClass : {
      std::string errorMessage = "Unsupported refinement functor type MultipleClass "
          "for density estimation!";
      throw new application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::Classification : {
      std::string errorMessage = "Unsupported refinement functor type Classification "
          "for density estimation!";
      throw new application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::GridPointBased : {
      std::string errorMessage = "Unsupported refinement functor type GridPointBased "
          "for density estimation!";
      throw new application_exception(errorMessage.c_str());
    }
  }
  return nullptr;
}

bool ModelFittingDensityEstimation::refine() {
  if (grid != nullptr && this->isRefinable()) {
    if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
      // create refinement functor
      RefinementFunctor *func = getRefinementFunctor();
      if (func != nullptr) {
        // refine grid
        auto oldNoPoints = grid->getSize();
        std::cout << "Old number points " << oldNoPoints << std::endl;
        grid->getGenerator().refine(*func);
        auto newNoPoints = grid->getSize();
        std::cout << "New number points " << newNoPoints << std::endl;
        if (newNoPoints != oldNoPoints) {
          // TODO(roehner) enable coarsening
          std::list<size_t> deletedGridPoints {};
          this->refine(newNoPoints, &deletedGridPoints);
          refinementsPerformed++;
          return true;
        } else {
          return false;
        }
      } else {
        throw application_exception(
            "ModelFittingDensityEstimation: No refinement functor could be created!");
      }
    } else {
      return false;
    }

  } else {
    throw application_exception(
        "ModelFittingDensityEstimation: Can't refine before initial grid is created");
  }
}

}  // namespace datadriven
}  // namespace sgpp
