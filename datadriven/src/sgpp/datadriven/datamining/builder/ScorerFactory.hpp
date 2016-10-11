/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ScorerFactory.hpp
 *
 * Created on: Oct 11, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/RandomShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/SequentialShufflingFunctor.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

class ScorerFactory {
 public:
  ScorerFactory(){};
  virtual ~ScorerFactory(){};

  virtual Scorer* buildScorer(const DataMiningConfigParser& parser) = 0;

 protected:
  std::unique_ptr<Metric> buildMetric(ScorerMetric config) {
    if (config == ScorerMetric::MSE) {
      return std::make_unique<MSE>();
    } else {
      return nullptr;
    }
  }

  std::unique_ptr<ShufflingFunctor> buildShuffling(ScorerShufflingType config) {
    if (config == ScorerShufflingType::random) {
      return std::make_unique<RandomShufflingFunctor>();
    } else if (config == ScorerShufflingType::sequential) {
      return std::make_unique<SequentialShufflingFunctor>();
    } else {
      return nullptr;
    }
  }
};

} /* namespace datadriven */
} /* namespace sgpp */
