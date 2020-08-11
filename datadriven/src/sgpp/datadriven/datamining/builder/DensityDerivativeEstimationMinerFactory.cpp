// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/DensityDerivativeEstimationMinerFactory.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
// #include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationCG.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityDerivativeEstimationOnOffParallel.hpp>
/*
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/DensityEstimationFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>
*/
// #include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDummy.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

ModelFittingBase *DensityDerivativeEstimationMinerFactory::createFitter(
    const DataMiningConfigParser &parser) const {
  FitterConfigurationDensityEstimation config{};
  config.readParams(parser);
  #ifdef USE_SCALAPACK
    if (parser.hasParallelConfig()) {
      return new ModelFittingDensityDerivativeEstimationOnOffParallel(config);
    }
  #endif
  /*
  if (config.getGridConfig().generalType_ == base::GeneralGridType::ComponentGrid) {
    return new ModelFittingDensityEstimationCombi(config);
  }
  */
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

/*
HyperparameterOptimizer *DensityDerivativeEstimationMinerFactory::buildHPO(
    const std::string &path) const {
  DataMiningConfigParser parser(path);
  if (parser.getHPOMethod("bayesian") == "harmonica") {
    return new HarmonicaHyperparameterOptimizer(buildMiner(path),
                                                new DensityEstimationFitterFactory(parser), parser);
  } else {
    return new BoHyperparameterOptimizer(buildMiner(path),
                                         new DensityEstimationFitterFactory(parser), parser);
  }
}
*/

/*
FitterFactory *DensityDerivativeEstimationMinerFactory::createFitterFactory(
    const DataMiningConfigParser &parser) const {
  return new DensityEstimationFitterFactory(parser);
}
*/

Visualizer *DensityDerivativeEstimationMinerFactory::createVisualizer(
    const DataMiningConfigParser &parser) const {
  VisualizerConfiguration config;
  config.readParams(parser);

  // TODO(spc90): implement visualization for this model
  return new VisualizerDummy();
}

} /* namespace datadriven */
} /* namespace sgpp */
