// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/DensityDerivativeRatioEstimationMinerFactory.hpp>

#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeRatioEstimation.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

ModelFittingBase* DensityDerivativeRatioEstimationMinerFactory::createFitter(
    const DataMiningConfigParser& parser) const {
  FitterConfigurationLeastSquares config{};
  config.readParams(parser);
  // Sanity check: will throw before creating anything if existence conditions are not met
  sanityCheck(config);
  return new ModelFittingDensityDerivativeRatioEstimation(config);
}

void DensityDerivativeRatioEstimationMinerFactory::sanityCheck(
    const FitterConfigurationLeastSquares& config) const {
  // Sanity check: SGDDerivRE needs Bspline basis of order at least 3
  if (config.getGridConfig().maxDegree_ < 3 &&
      !(config.getGridConfig().type_ == base::GridType::Bspline ||
        config.getGridConfig().type_ == base::GridType::BsplineBoundary))
    throw base::algorithm_exception(
        "DensityDerivativeRatioEstimationMinerFactory: Method requires Bspline basis functions of "
        "order at least 3!");
}

/*
FitterFactory* DensityDerivativeRatioEstimationMinerFactory::createFitterFactory(
    const DataMiningConfigParser& parser) const {
  return new LeastSquaresRegressionFitterFactory(parser);
}
*/

Visualizer* DensityDerivativeRatioEstimationMinerFactory::createVisualizer(
    const DataMiningConfigParser& parser) const {
  VisualizerConfiguration config;
  config.readParams(parser);

  return new VisualizerDensityEstimation(config);
}

} /* namespace datadriven */
} /* namespace sgpp */
