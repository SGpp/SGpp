// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/NegativeLogLikelihood.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Accuracy.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ResidualScore.hpp>

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
    // (Sebastian Kreisel) This case should never occur, because
    // ScorerConfiguration sets up a default value: ScorerMetricType:accuracy
    // Previously we returned a nullptr in this else-case but this leads
    // to segfaults or undefined behavior down the line, so I removed it
    // TODO(Sebastian) It would be best to throw an error or exception here
    return new Accuracy{};
  }
}

Scorer* ScorerFactory::buildScorer(const DataMiningConfigParser& parser) {
  ScorerConfiguration config;
  parser.getScorerConfig(config, config);
  auto metric = buildMetric(config.metric_);
  return new Scorer(metric);
}

Metric* ScorerFactory::buildRegularizationMetric(RegularizationMetricType config) const {
  if (config == RegularizationMetricType::mse) {
    return new MSE{};
  } else if (config == RegularizationMetricType::nll) {
    return new NegativeLogLikelihood{};
  } else if (config == RegularizationMetricType::accuracy) {
    return new Accuracy{};
  } else if (config == RegularizationMetricType::residual) {
    return new ResidualScore{};
  } else {
    // (Sebastian Kreisel) This case should never occur, because
    // ScorerConfiguration sets up a default value: RegularizationMetricType::Residual
    // Previously we returned a nullptr in this else-case but this leads
    // to segfaults or undefined behavior down the line, so I removed it
    // TODO(Sebastian) It would be best to throw an error or exception here
    return new ResidualScore{};
  }
}

Scorer* ScorerFactory::buildRegularizationScorer(const RegularizationConfiguration& config) {
  auto metric = buildRegularizationMetric(config.regularizationMetric_);
  return new Scorer(metric);
}

} /* namespace datadriven */
} /* namespace sgpp */
