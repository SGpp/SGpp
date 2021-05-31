// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>

#include <list>
#include <string>
#include <vector>

using sgpp::base::CoarseningFunctor;
using sgpp::base::CoarseningFunctorType;
using sgpp::base::RefinementFunctor;
using sgpp::base::RefinementFunctorType;
using sgpp::base::SurplusCoarseningFunctor;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::SurplusVolumeCoarseningFunctor;
using sgpp::base::SurplusVolumeRefinementFunctor;

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimation::ModelFittingDensityEstimation()
    : refinementsPerformed(0), initialGridSize(0) {}

std::unique_ptr<RefinementFunctor> ModelFittingDensityEstimation::getRefinementFunctor() {
  sgpp::base::AdaptivityConfiguration &refinementConfig = this->config->getRefinementConfig();
  switch (refinementConfig.refinementFunctorType_) {
    case RefinementFunctorType::Surplus: {
      return std::make_unique<SurplusRefinementFunctor>(
          alpha, config->getRefinementConfig().numRefinementPoints_,
          config->getRefinementConfig().refinementThreshold_);
    }
    case RefinementFunctorType::SurplusVolume: {
      return std::make_unique<SurplusVolumeRefinementFunctor>(
          alpha, config->getRefinementConfig().numRefinementPoints_,
          config->getRefinementConfig().refinementThreshold_);
    }
    case RefinementFunctorType::DataBased: {
      std::string errorMessage =
          "Unsupported refinement functor type DataBased "
          "for density estimation!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::ZeroCrossing: {
      std::string errorMessage =
          "Unsupported refinement functor type ZeroCrossing "
          "for density estimation!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::MultipleClass: {
      std::string errorMessage =
          "Unsupported refinement functor type MultipleClass "
          "for density estimation!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::Classification: {
      std::string errorMessage =
          "Unsupported refinement functor type Classification "
          "for density estimation!";
      throw application_exception(errorMessage.c_str());
    }
    case RefinementFunctorType::GridPointBased: {
      std::string errorMessage =
          "Unsupported refinement functor type GridPointBased "
          "for density estimation!";
      throw application_exception(errorMessage.c_str());
    }
  }
  return nullptr;
}

std::unique_ptr<CoarseningFunctor> ModelFittingDensityEstimation::getCoarseningFunctor() {
  sgpp::base::AdaptivityConfiguration &adaptivityConfig = this->config->getRefinementConfig();
  switch (adaptivityConfig.coarseningFunctorType_) {
    case CoarseningFunctorType::Surplus: {
      return std::make_unique<SurplusCoarseningFunctor>(
          alpha, config->getRefinementConfig().numRefinementPoints_,
          config->getRefinementConfig().coarseningThreshold_);
    }
    case CoarseningFunctorType::SurplusVolume: {
      return std::make_unique<SurplusVolumeCoarseningFunctor>(
          alpha, config->getRefinementConfig().numRefinementPoints_,
          config->getRefinementConfig().coarseningThreshold_);
    }
    case CoarseningFunctorType::Classification: {
      std::string errorMessage =
          "Unsupported refinement functor type Classification "
          "for density estimation!";
      throw application_exception(errorMessage.c_str());
    }
  }
  return nullptr;
}

bool ModelFittingDensityEstimation::adapt() {
  if (grid != nullptr && this->isRefinable()) {
    if (this->initialGridSize == 0) {
      this->initialGridSize = grid->getSize();
    }

    if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
      // create refinement and coarsening functors
      std::unique_ptr<RefinementFunctor> refinementFunc = getRefinementFunctor();
      auto oldNoPoints = grid->getSize();

      std::unique_ptr<CoarseningFunctor> coarseningFunc = getCoarseningFunctor();
      std::vector<size_t> deletedGridPoints;

      // do coarsening before refinement to prevent refined grid points from being coarsened
      // immediately

      if (coarseningFunc) {
        // coarsen grid
        if (config->getRefinementConfig().coarsenInitialPoints_) {
          grid->getGenerator().coarsenNFirstOnly(*coarseningFunc, grid->getSize(),
                                                 &deletedGridPoints, 0);
        } else {
          grid->getGenerator().coarsenNFirstOnly(*coarseningFunc, grid->getSize(),
                                                 &deletedGridPoints, this->initialGridSize);
        }
      } else {
        throw application_exception(
            "ModelFittingDensityEstimation: No coarsening functor could be "
            "created!");
      }

      if (refinementFunc) {
        // refine grid
        std::cout << "Old number points " << oldNoPoints << std::endl;
        GeometryConfiguration geometryConfig = config->getGeometryConfig();
        if (!geometryConfig.stencils_.empty()) {
          GridFactory gridFactory;
          grid->getGenerator().refineInter(*refinementFunc,
                                           gridFactory.getInteractions(geometryConfig));
        } else {
          grid->getGenerator().refine(*refinementFunc);
        }
      } else {
        throw application_exception(
            "ModelFittingDensityEstimation: No refinement functor could be "
            "created!");
      }

      auto newNoPoints = grid->getSize();
      std::cout << "New number points " << newNoPoints << std::endl;
      if (newNoPoints != oldNoPoints || !deletedGridPoints.empty()) {
        this->adapt(newNoPoints, deletedGridPoints);
        refinementsPerformed++;
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else {
    throw application_exception(
        "ModelFittingDensityEstimation: Can't refine before initial grid is "
        "created");
  }
}

}  // namespace datadriven
}  // namespace sgpp
