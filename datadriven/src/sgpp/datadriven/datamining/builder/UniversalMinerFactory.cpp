/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * UniversalMinerFactory.cpp
 *
 * Created on: Mar 12, 2018
 *     Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/builder/UniversalMinerFactory.hpp>

#include <string>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/CrossValidationScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/SplittingScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/LeastSquaresRegressionFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/DensityEstimationFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>

namespace sgpp {
namespace datadriven {

SparseGridMiner *UniversalMinerFactory::buildMiner(const std::string &path) const {
  DataMiningConfigParser parser(path);

  return new SparseGridMiner(createDataSource(parser), createFitter(parser), createScorer(parser));
}

HyperparameterOptimizer *UniversalMinerFactory::buildHPO(const std::string &path) const {
  DataMiningConfigParser parser(path);
  FitterType fType = FitterType::RegressionLeastSquares;
  parser.getFitterConfigType(fType, fType);

  if (fType == FitterType::DensityEstimation) {
    return new HyperparameterOptimizer(createDataSource(parser),
                                       new DensityEstimationFitterFactory(parser),
                                       parser);
  } else {
    return new HyperparameterOptimizer(createDataSource(parser),
                                       new LeastSquaresRegressionFitterFactory(parser),
                                       parser);
  }
}

DataSource *UniversalMinerFactory::createDataSource(
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

ModelFittingBase *UniversalMinerFactory::createFitter(
    const DataMiningConfigParser &parser) const {
  ModelFittingBase *model;

  FitterType fType = FitterType::RegressionLeastSquares;
  parser.getFitterConfigType(fType, fType);

  if (fType == FitterType::DensityEstimation) {
    FitterConfigurationDensityEstimation config{};
    config.readParams(parser);
    model = new ModelFittingDensityEstimation(config);
  } else {
    FitterConfigurationLeastSquares config{};
    config.readParams(parser);
    model = new ModelFittingLeastSquares(config);
  }
  return model;
}

Scorer *UniversalMinerFactory::createScorer(
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
