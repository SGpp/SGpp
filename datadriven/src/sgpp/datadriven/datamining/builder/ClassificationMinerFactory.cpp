/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ClassificationMinerFactory.cpp
 *
 *  Created on: Jul 2, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/CrossValidationScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/SplittingScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

SparseGridMiner* ClassificationMinerFactory::buildMiner(const std::string& path) const {
  DataMiningConfigParser parser(path);

  return new SparseGridMiner(createDataSource(parser), createFitter(parser), createScorer(parser));
}

DataSource* ClassificationMinerFactory::createDataSource(
    const DataMiningConfigParser& parser) const {
  DataSourceConfig config;

  bool hasSource = parser.getDataSourceConfig(config, config);

  if (hasSource && config.filePath.compare("") != 0) {
    DataSourceBuilder builder;
    return builder.fromConfig(config);
  } else {
    throw base::data_exception("No file name provided for datasource.");
  }
}

ModelFittingBase* ClassificationMinerFactory::createFitter(
    const DataMiningConfigParser& parser) const {
  FitterConfigurationClassification config{};
  config.readParams(parser);
  return new ModelFittingClassification(config);
}

Scorer* ClassificationMinerFactory::createScorer(
    const DataMiningConfigParser& parser) const {
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
