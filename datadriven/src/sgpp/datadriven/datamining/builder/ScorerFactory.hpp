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
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Abstract factory to build all kinds of scorers based on a given configuration.
 */
class ScorerFactory {
 public:
  /**
   * Default constructor
   */
  ScorerFactory() = default;

  /**
   * Virtual destructor
   */
  virtual ~ScorerFactory() = default;

  /**
   * Assemble a #sgpp::datadriven::Scorer object based on the configuration
   * @param parser Instance of #sgpp::datadriven::DataMiningConfigParser that reads the required
   * data from the config file.
   * @return Fully configured instance of a  #sgpp::datadriven::Scorer object.
   */
  virtual Scorer* buildScorer(const DataMiningConfigParser& parser) const = 0;

 protected:
  /**
   * Build a #sgpp::datadriven::Metric object based on the given metric type enum value.
   * @param config #sgpp::datadriven::ScorerMetricType describing which #sgpp::datadriven::Metric to
   * generate.
   * @return  Fully configured instance of a  #sgpp::datadriven::Metric object.
   */
  Metric* buildMetric(ScorerMetricType config) const;

  /**
   * Build a #sgpp::datadriven::ShufflingFunctor object based on the the given shuffling type enum
   * value.
   * @param config #sgpp::datadriven::ScorerShufflingType describing which
   * #sgpp::datadriven::ShufflingFunctor to generate.
   * @return Fully configured instance of a  #sgpp::datadriven::ShufflingFunctor object.
   */
  ShufflingFunctor* buildShuffling(ScorerShufflingType config) const;
};
} /* namespace datadriven */
} /* namespace sgpp */
