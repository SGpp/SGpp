// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/algorithm/DBMatObjectStore.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Class to add combigrid functionality to density derivative estimators.
 */

class ModelFittingDensityDerivativeEstimationCombi : public ModelFittingDensityEstimationCombi {
 public:
  /**
   * Default constructor
   */
  ModelFittingDensityDerivativeEstimationCombi();

  /**
   * Constructor from a configuration object generated by the datamining
   * pipeline
   * @param config configuration object generated by the datamining pipeline
   */
  explicit ModelFittingDensityDerivativeEstimationCombi(
      const FitterConfigurationDensityEstimation& config);

  /**
   * @brief Construct from a configuration object generated from the datamining
   * pipeline and a
   * object store to obtain and store already decomposed offline objects.
   *
   * @param config   Configuration object generated by the datamining pipeline.
   * @param objectStore  Offline object store.
   */
  explicit ModelFittingDensityDerivativeEstimationCombi(
      const FitterConfigurationDensityEstimation& config,
      std::shared_ptr<DBMatObjectStore> objectStore);

 protected:
  /**
   * Creates a density estimation model that fits the model settings.
   * @param densityEstimationConfig configuration for the density estimation
   * @return a new density estimation model
   */
  std::unique_ptr<ModelFittingDensityEstimation> createNewModel(
      sgpp::datadriven::FitterConfigurationDensityEstimation& densityEstimationConfig) override;
};

}  // namespace datadriven
}  // namespace sgpp