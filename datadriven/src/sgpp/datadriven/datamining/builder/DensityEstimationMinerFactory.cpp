// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/DensityEstimationMinerFactory.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/DensityEstimationFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

ModelFittingBase *DensityEstimationMinerFactory::createFitter(
    const DataMiningConfigParser &parser) const {
  FitterConfigurationDensityEstimation config{};
  config.readParams(parser);
#ifdef USE_SCALAPACK
  if (parser.hasParallelConfig()) {
    return new ModelFittingDensityEstimationOnOffParallel(config);
  }
#endif
  if (config.getGridConfig().generalType_ == base::GeneralGridType::ComponentGrid) {
    return new ModelFittingDensityEstimationCombi(config);
  }
  switch (config.getDensityEstimationConfig().type_) {
    case (DensityEstimationType::CG):
      std::cout << "\nCG\n";
      return new ModelFittingDensityEstimationCG(config);
    case (DensityEstimationType::Decomposition):
      std::cout << "\nDECOMPOSITION\n";
      return new ModelFittingDensityEstimationOnOff(config);
  }

  throw base::application_exception("Unknown density estimation type");
}

HyperparameterOptimizer *DensityEstimationMinerFactory::buildHPO(const std::string &path) const {
  DataMiningConfigParser parser(path);
  if (parser.getHPOMethod("bayesian") == "harmonica") {
    return new HarmonicaHyperparameterOptimizer(buildMiner(path),
                                                new DensityEstimationFitterFactory(parser), parser);
  } else {
    return new BoHyperparameterOptimizer(buildMiner(path),
                                         new DensityEstimationFitterFactory(parser), parser);
  }
}
FitterFactory *DensityEstimationMinerFactory::createFitterFactory(
    const DataMiningConfigParser &parser) const {
  return new DensityEstimationFitterFactory(parser);
}

Visualizer *DensityEstimationMinerFactory::createVisualizer(
    const DataMiningConfigParser &parser) const {
  VisualizerConfiguration config;
  config.readParams(parser);

  return new VisualizerDensityEstimation(config);
}

} /* namespace datadriven */
} /* namespace sgpp */
