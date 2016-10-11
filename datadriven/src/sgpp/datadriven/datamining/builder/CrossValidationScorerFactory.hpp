/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossValidationScorerFactory.hpp
 *
 * Created on: Oct 11, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include "ScorerFactory.hpp"

#include <sgpp/datadriven/datamining/modules/scoring/CrossValidation.hpp>

namespace sgpp {
namespace datadriven {

class CrossValidationScorerFactory : public ScorerFactory {
 public:
  CrossValidationScorerFactory() : ScorerFactory(){};
  virtual ~CrossValidationScorerFactory(){};

  virtual Scorer* buildScorer(const DataMiningConfigParser& parser) {
    CrossValidationConfiguration config;
    parser.getScorerCrossValidationConfig(config, config);

    auto metric = buildMetric(config.metric);
    auto shuffling = buildShuffling(config.shuffling);

    return new CrossValidation(metric.release(), shuffling.release(), config.randomSeed,
                               config.folds);
  };
};

} /* namespace datadriven */
} /* namespace sgpp */
