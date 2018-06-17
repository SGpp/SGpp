/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * HPOScorerFactory.cpp
 *
 *  Created on:	10.12.2017
 *      Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/modules/hpo/HPOScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HPOScorer.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>

namespace sgpp {
namespace datadriven {

Scorer *HPOScorerFactory::buildScorer(const DataMiningConfigParser &parser) const {
  TestingConfiguration config;
  parser.getScorerTestingConfig(config, config);
  auto metric = buildMetric(config.metric);
  auto shuffling = buildShuffling(config.shuffling);

  DataSourceConfig dsConfig;

  bool hasSource = parser.getScorerTestset(dsConfig, dsConfig);

  if ((hasSource && dsConfig.filePath.compare("") != 0)) {

    DataSourceBuilder builder;
    std::unique_ptr<DataSource> dataSource;
    dataSource.reset(builder.fromConfig(dsConfig));

    return new HPOScorer(metric,
                         shuffling,
                         config.randomSeed,
                         config.testingPortion,
                         dataSource->getNextSamples());
  }

  return new HPOScorer(metric, shuffling, config.randomSeed, config.testingPortion, nullptr);
}
} /* namespace datadriven */
} /* namespace sgpp */
