// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/DensityDerivativeEstimationMinerFactory.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationCombi.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationOnOffParallel.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/DensityDerivativeEstimationFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

ModelFittingBase *DensityDerivativeEstimationMinerFactory::createFitter(
    const DataMiningConfigParser &parser) const {
  FitterConfigurationDensityEstimation config{};
  config.readParams(parser);
  // Sanity check: will throw before creating anything if existence conditions are not met
  sanityCheck(config);
#ifdef USE_SCALAPACK
  if (parser.hasParallelConfig()) {
    return new ModelFittingDensityDerivativeEstimationOnOffParallel(config);
  }
#endif
  if (config.getGridConfig().generalType_ == base::GeneralGridType::ComponentGrid) {
    return new ModelFittingDensityDerivativeEstimationCombi(config);
  }
  switch (config.getDensityEstimationConfig().type_) {
    case (DensityEstimationType::CG):
      std::cout << "\nCG\n";
      return new ModelFittingDensityDerivativeEstimationCG(config);
    case (DensityEstimationType::Decomposition):
      std::cout << "\nDECOMPOSITION\n";
      return new ModelFittingDensityDerivativeEstimationOnOff(config);
  }

  throw base::application_exception("Unknown density estimation type");
}

void DensityDerivativeEstimationMinerFactory::sanityCheck(
    const FitterConfigurationDensityEstimation &config) const {
  // Sanity check: SGDDerivRE needs Bspline basis of order at least 3
  if (config.getGridConfig().maxDegree_ < 3 &&
      !(config.getGridConfig().type_ == base::GridType::Bspline ||
        config.getGridConfig().type_ == base::GridType::BsplineBoundary ||
        config.getGridConfig().type_ == base::GridType::ModBspline ||
        config.getGridConfig().type_ == base::GridType::NakBsplineBoundary ||
        config.getGridConfig().type_ == base::GridType::NakBsplineExtended ||
        config.getGridConfig().type_ == base::GridType::ModNakBspline))
    throw base::algorithm_exception(
        "DensityDerivativeEstimationMinerFactory: Method requires Bspline basis functions of "
        "degree at least 3!");
}

HyperparameterOptimizer *DensityDerivativeEstimationMinerFactory::buildHPO(
    const std::string &path) const {
  DataMiningConfigParser parser(path);
  if (parser.getHPOMethod("bayesian") == "harmonica") {
    return new HarmonicaHyperparameterOptimizer(
        buildMiner(path), new DensityDerivativeEstimationFitterFactory(parser), parser);
  } else {
    return new BoHyperparameterOptimizer(
        buildMiner(path), new DensityDerivativeEstimationFitterFactory(parser), parser);
  }
}

FitterFactory *DensityDerivativeEstimationMinerFactory::createFitterFactory(
    const DataMiningConfigParser &parser) const {
  return new DensityDerivativeEstimationFitterFactory(parser);
}

Visualizer *DensityDerivativeEstimationMinerFactory::createVisualizer(
    const DataMiningConfigParser &parser) const {
  VisualizerConfiguration config;
  config.readParams(parser);

  return new VisualizerDensityEstimation(config);
}

} /* namespace datadriven */
} /* namespace sgpp */
