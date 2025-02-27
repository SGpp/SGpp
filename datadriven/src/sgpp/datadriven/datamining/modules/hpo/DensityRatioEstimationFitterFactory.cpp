// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/DensityRatioEstimationFitterFactory.hpp>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityRatioEstimation.hpp>

namespace sgpp {
namespace datadriven {

DensityRatioEstimationFitterFactory::DensityRatioEstimationFitterFactory(
    const DataMiningConfigParser &parser)
    : baseConfig() {
  baseConfig.readParams(parser);

  parser.getHyperparameters(conpar, dispar, catpar, basisFunctions);
}

ModelFittingBase *DensityRatioEstimationFitterFactory::buildFitter() {
  // build config
  auto *config = new FitterConfigurationLeastSquares(baseConfig);

  if (dispar.count("level")) {
    config->getGridConfig().level_ = dispar["level"].getValue();
  }
  if (catpar.count("basisFunction")) {
    config->getGridConfig().type_ = basisFunctions[catpar["basisFunction"].getValue()];
  }
  if (dispar.count("noPoints")) {
    config->getRefinementConfig().numRefinementPoints_ =
        static_cast<size_t>(dispar["noPoints"].getValue());
  }
  if (conpar.count("threshold")) {
    config->getRefinementConfig().refinementThreshold_ = conpar["threshold"].getValue();
  }
  if (conpar.count("lambda")) {
    config->getRegularizationConfig().lambda_ = conpar["lambda"].getValue();
  }

  return new ModelFittingDensityRatioEstimation(*config);
}
} /* namespace datadriven */
} /* namespace sgpp */
