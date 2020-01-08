// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>

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

ModelFittingDensityEstimation::ModelFittingDensityEstimation() : refinementsPerformed{0} {}

std::unique_ptr<RefinementFunctor> ModelFittingDensityEstimation::getRefinementFunctor() {
  sgpp::base::AdaptivityConfiguration &refinementConfig = this->config->getRefinementConfig();
  switch (refinementConfig.refinementFunctorType) {
    case RefinementFunctorType::Surplus: {
      return std::make_unique<SurplusRefinementFunctor>(
          alpha, config->getRefinementConfig().noPoints_,
          config->getRefinementConfig().refinementThreshold_);
    }
    case RefinementFunctorType::SurplusVolume: {
      return std::make_unique<SurplusVolumeRefinementFunctor>(
          alpha, config->getRefinementConfig().noPoints_,
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
  switch (adaptivityConfig.coarseningFunctorType) {
    case CoarseningFunctorType::Surplus: {
      return std::make_unique<SurplusCoarseningFunctor>(
          alpha, config->getRefinementConfig().noPoints_,
          config->getRefinementConfig().coarseningThreshold_);
    }
    case CoarseningFunctorType::SurplusVolume: {
      return std::make_unique<SurplusVolumeCoarseningFunctor>(
          alpha, config->getRefinementConfig().noPoints_,
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
    if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
      // create refinement functor
      std::unique_ptr<RefinementFunctor> refinementFunc = getRefinementFunctor();
      auto oldNoPoints = grid->getSize();
      if (refinementFunc) {
        // refine grid
        std::cout << "Old number points " << oldNoPoints << std::endl;
        GeometryConfiguration geoConf = config->getGeometryConfig();
        if (!geoConf.stencils.empty()) {
          GridFactory gridFactory;
          grid->getGenerator().refineInter(*refinementFunc, gridFactory.getInteractions(geoConf));
        } else {
          grid->getGenerator().refine(*refinementFunc);
        }
      } else {
        throw application_exception(
            "ModelFittingDensityEstimation: No refinement functor could be created!");
      }

      std::unique_ptr<CoarseningFunctor> coarseningFunc = getCoarseningFunctor();
      std::list<size_t> deletedGridPoints;

      // coarsen grid
      if (coarseningFunc) {
        // TODO(jan schopohl) coarsen n first only?
        std::vector<size_t> removedSeq;
        grid->getGenerator().coarsen(*coarseningFunc, alpha, &removedSeq);
        std::copy(removedSeq.begin(), removedSeq.end(), std::back_inserter(deletedGridPoints));
      } else {
        // TODO (jan schopohl) throw exception or just don't do coarsening/refine?
        throw application_exception(
            "ModelFittingDensityEstimation: No coarsening functor could be created!");
      }

      auto newNoPoints = grid->getSize();
      std::cout << "New number points " << newNoPoints << std::endl;
      if (newNoPoints != oldNoPoints) {
        this->adapt(newNoPoints, &deletedGridPoints);
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
        "ModelFittingDensityEstimation: Can't refine before initial grid is created");
  }
}  // namespace datadriven

}  // namespace datadriven
}  // namespace sgpp
