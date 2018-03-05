/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SplittingScorerFactory.cpp
 *
 *  Created on: 25.01.2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/builder/SplittingScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/SplittingScorer.hpp>

namespace sgpp {
namespace datadriven {
Scorer* SplittingScorerFactory::buildScorer(const DataMiningConfigParser& parser) const {
  TestingConfiguration config;
  parser.getScorerTestingConfig(config, config);

  auto metric = buildMetric(config.metric);
  auto shuffling = buildShuffling(config.shuffling);
  return new SplittingScorer(metric, shuffling, config.randomSeed, config.testingPortion);
}
} /* namespace datadriven */
} /* namespace sgpp */
