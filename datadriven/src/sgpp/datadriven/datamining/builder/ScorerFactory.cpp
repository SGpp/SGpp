/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ScorerFactory.cpp
 *
 *  Created on: 25.01.2017
 *      Author: Michael Lettrich
 */
#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/NegativeLogLikelihood.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Accuracy.hpp>

namespace sgpp {
namespace datadriven {
Metric* ScorerFactory::buildMetric(ScorerMetricType config) const {
  if (config == ScorerMetricType::mse) {
    return new MSE{};
  } else if (config == ScorerMetricType::nll) {
    return new NegativeLogLikelihood{};
  } else if (config == ScorerMetricType::accuracy) {
    return new Accuracy{};
  } else {
    return nullptr;
  }
}

Scorer* ScorerFactory::buildScorer(const DataMiningConfigParser& parser) {
  ScorerConfiguration config;
  parser.getScorerConfig(config, config);
  auto metric = buildMetric(config.metric);
  return new Scorer(metric);
}


} /* namespace datadriven */
} /* namespace sgpp */
