// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Factory to build the scorer
 */
class ScorerFactory {
 public:
  /**
   * Default constructor
   */
  ScorerFactory() = default;

  /**
   * Destructor
   */
  ~ScorerFactory() = default;

  /**
   * Assemble a #sgpp::datadriven::Scorer object based on the configuration
   * @param parser Instance of #sgpp::datadriven::DataMiningConfigParser that reads the required
   * data from the config file.
   * @return Fully configured instance of a  #sgpp::datadriven::Scorer object.
   */
  Scorer* buildScorer(const DataMiningConfigParser& parser);

  /**
   * Assemble a #sgpp::datadriven::Scorer for regularization optimization based on the configuration
   * @param config regularization config
   * @return Fully configured instance of a  #sgpp::datadriven::Scorer object.
   */
  Scorer* buildRegularizationScorer(const RegularizationConfiguration& config);

 protected:
  /**
   * Build a #sgpp::datadriven::Metric object based on the given metric type enum value.
   * @param config #sgpp::datadriven::ScorerMetricType describing which #sgpp::datadriven::Metric to
   * generate.
   * @return  Fully configured instance of a  #sgpp::datadriven::Metric object.
   */
  Metric* buildMetric(ScorerMetricType config) const;

  /**
 * Build a #sgpp::datadriven::Metric object based on the given metric type enum value.
 * @param config #sgpp::datadriven::RegularizationMetricType describing which
 * #sgpp::datadriven::Metric to
 * generate.
 * @return  Fully configured instance of a  #sgpp::datadriven::Metric object.
 */
  Metric* buildRegularizationMetric(RegularizationMetricType config) const;
};
} /* namespace datadriven */
} /* namespace sgpp */
