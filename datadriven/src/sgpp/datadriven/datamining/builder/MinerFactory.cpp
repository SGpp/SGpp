// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/MinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerCrossValidation.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMinerSplitting.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

SparseGridMiner* MinerFactory::buildMiner(const std::string& path) const {
  DataMiningConfigParser parser(path);

  if (parser.hasFitterConfigCrossValidation()) {
     return new SparseGridMinerCrossValidation(createDataSourceCrossValidation(parser),
                                               createFitter(parser), createScorer(parser),
                                               createVisualizer(parser));
  } else {
     return new SparseGridMinerSplitting(createDataSourceSplitting(parser), createFitter(parser),
                                         createScorer(parser),
                                         createVisualizer(parser));
  }
}

sgpp::datadriven::HyperparameterOptimizer* MinerFactory::buildHPO(const std::string& path) const {
  DataMiningConfigParser parser(path);
  if (parser.getHPOMethod("bayesian") == "harmonica") {
    return new HarmonicaHyperparameterOptimizer(buildMiner(path), createFitterFactory(parser),
                                                parser);
  } else {
    return new BoHyperparameterOptimizer(buildMiner(path), createFitterFactory(parser), parser);
  }
}

DataSourceSplitting* MinerFactory::createDataSourceSplitting(
    const DataMiningConfigParser& parser) const {
  DataSourceConfig config{};

  bool hasSource = parser.getDataSourceConfig(config, config);

  if (hasSource && config.filePath_.compare("") != 0) {
    DataSourceBuilder builder;
    return builder.splittingFromConfig(config);
  } else {
    throw base::data_exception("No file name provided for datasource.");
  }
}

DataSourceCrossValidation* MinerFactory::createDataSourceCrossValidation(
    const DataMiningConfigParser& parser) const {
  DataSourceConfig config{};

  CrossvalidationConfiguration crossValidationconfig{};
  parser.getFitterCrossvalidationConfig(crossValidationconfig, crossValidationconfig);

  bool hasSource = parser.getDataSourceConfig(config, config);

  if (hasSource && config.filePath_.compare("") != 0) {
    DataSourceBuilder builder;
    return builder.crossValidationFromConfig(config, crossValidationconfig);
  } else {
    throw base::data_exception("No file name provided for datasource.");
  }
}

Scorer* MinerFactory::createScorer(const DataMiningConfigParser& parser) const {
  std::unique_ptr<ScorerFactory> factory = std::make_unique<ScorerFactory>();
  return factory->buildScorer(parser);
}


} /* namespace datadriven */
} /* namespace sgpp */
