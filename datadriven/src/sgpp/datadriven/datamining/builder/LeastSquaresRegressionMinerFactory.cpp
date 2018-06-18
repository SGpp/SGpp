/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * LeastSquaresRegressionFactory.cpp
 *
 * Created on: Oct 10, 2016
 *     Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/builder/LeastSquaresRegressionMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/CrossValidationScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/SplittingScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/LeastSquaresRegressionFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

SparseGridMiner *LeastSquaresRegressionMinerFactory::buildMiner(const std::string &path) const {
  DataMiningConfigParser parser(path);

  return new SparseGridMiner(createDataSource(parser), createFitter(parser), createScorer(parser));
}

HyperparameterOptimizer *LeastSquaresRegressionMinerFactory::buildHPO
                                                             (const std::string &path) const {
  DataMiningConfigParser parser(path);
  return new HyperparameterOptimizer(createDataSource(parser),
                                     new LeastSquaresRegressionFitterFactory(parser),
                                     parser);
}

DataSource *LeastSquaresRegressionMinerFactory::createDataSource(
    const DataMiningConfigParser &parser) const {
  DataSourceConfig config;

  bool hasSource = parser.getDataSourceConfig(config, config);

  if (hasSource && config.filePath.compare("") != 0) {
    DataSourceBuilder builder;
    return builder.fromConfig(config);
  } else {
    throw base::data_exception("No file name provided for datasource.");
  }
}

ModelFittingBase *LeastSquaresRegressionMinerFactory::createFitter(
    const DataMiningConfigParser &parser) const {
  FitterConfigurationLeastSquares config{};
  config.readParams(parser);
  return new ModelFittingLeastSquares(config);
}

Scorer *LeastSquaresRegressionMinerFactory::createScorer(
    const DataMiningConfigParser &parser) const {
  std::unique_ptr<ScorerFactory> factory;

  if (parser.hasScorerConfigCrossValidation()) {
    factory = std::make_unique<CrossValidationScorerFactory>();
  } else {
    factory = std::make_unique<SplittingScorerFactory>();
  }
  return factory->buildScorer(parser);
}
} /* namespace datadriven */
} /* namespace sgpp */
