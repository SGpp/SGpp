// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationCombi.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationOnOff.hpp>

using sgpp::base::application_exception;
using std::unique_ptr;
using std::vector;

namespace sgpp {
namespace datadriven {

ModelFittingDensityDerivativeEstimationCombi::ModelFittingDensityDerivativeEstimationCombi() {}

ModelFittingDensityDerivativeEstimationCombi::ModelFittingDensityDerivativeEstimationCombi(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimationCombi(config) {}

ModelFittingDensityDerivativeEstimationCombi::ModelFittingDensityDerivativeEstimationCombi(
    const FitterConfigurationDensityEstimation& config,
    std::shared_ptr<DBMatObjectStore> objectStore)
    : ModelFittingDensityEstimationCombi(config, objectStore) {}

std::unique_ptr<ModelFittingDensityEstimation>
ModelFittingDensityDerivativeEstimationCombi::createNewModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) {
  switch (densityEstimationConfig.getDensityEstimationConfig().type_) {
    case DensityEstimationType::CG: {
      return std::make_unique<ModelFittingDensityDerivativeEstimationCG>(densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
      if (this->hasObjectStore) {
        return std::make_unique<ModelFittingDensityDerivativeEstimationOnOff>(
            densityEstimationConfig, objectStore);
      } else {
        return std::make_unique<ModelFittingDensityDerivativeEstimationOnOff>(
            densityEstimationConfig);
      }
    }
  }

  throw base::application_exception("Unknown density estimation type");
}

}  // namespace datadriven
}  // namespace sgpp
