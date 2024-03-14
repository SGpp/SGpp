// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDifferenceEstimationCombi.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDifferenceEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDifferenceEstimationOnOff.hpp>

namespace sgpp {
namespace datadriven {

ModelFittingDensityDifferenceEstimationCombi::ModelFittingDensityDifferenceEstimationCombi() {}

ModelFittingDensityDifferenceEstimationCombi::ModelFittingDensityDifferenceEstimationCombi(
    const FitterConfigurationDensityEstimation& config)
    : ModelFittingDensityEstimationCombi(config) {}

ModelFittingDensityDifferenceEstimationCombi::ModelFittingDensityDifferenceEstimationCombi(
    const FitterConfigurationDensityEstimation& config,
    std::shared_ptr<DBMatObjectStore> objectStore)
    : ModelFittingDensityEstimationCombi(config, objectStore) {}

std::unique_ptr<ModelFittingDensityEstimation>
ModelFittingDensityDifferenceEstimationCombi::createNewModel(
    sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) {
  switch (densityEstimationConfig.getDensityEstimationConfig().type_) {
    case DensityEstimationType::CG: {
      return std::make_unique<ModelFittingDensityDifferenceEstimationCG>(densityEstimationConfig);
    }
    case DensityEstimationType::Decomposition: {
      return std::make_unique<ModelFittingDensityDifferenceEstimationOnOff>(
          densityEstimationConfig);
    }
  }

  throw base::application_exception("Unknown density estimation type");
}

}  // namespace datadriven
}  // namespace sgpp
