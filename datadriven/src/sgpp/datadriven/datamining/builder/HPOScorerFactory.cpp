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

#include <sgpp/datadriven/datamining/builder/HPOScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/HPOScorer.hpp>

namespace sgpp {
namespace datadriven {
Scorer* HPOScorerFactory::buildHPOScorer(const DataMiningConfigParser& parser, FitterFactory* fitterFactory) const {
  TestingConfiguration config;
  parser.getScorerTestingConfig(config, config);
  auto metric = buildMetric(config.metric);
  auto shuffling = buildShuffling(config.shuffling);
  return new HPOScorer(metric, shuffling, config.randomSeed, config.testingPortion, parser, fitterFactory);
}
Scorer* HPOScorerFactory::buildScorer(const DataMiningConfigParser& parser) const {
}
} /* namespace datadriven */
} /* namespace sgpp */
