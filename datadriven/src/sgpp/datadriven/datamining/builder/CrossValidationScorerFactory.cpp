/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossValidationScorerFactory.cpp
 *
 *  Created on: 25.01.2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/datamining/builder/CrossValidationScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/CrossValidation.hpp>

namespace sgpp {
namespace datadriven {
Scorer* CrossValidationScorerFactory::buildScorer(const DataMiningConfigParser& parser) const {
  CrossValidationConfiguration config;
  parser.getScorerCrossValidationConfig(config, config);

  auto metric = buildMetric(config.metric);
  auto shuffling = buildShuffling(config.shuffling);

  return new CrossValidation(metric, shuffling, config.randomSeed, config.folds);
}
} /* namespace datadriven */
} /* namespace sgpp */
