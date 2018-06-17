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

#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Concrete factory to build an instance of #sgpp::datadriven::CrossValidation
 */
class CrossValidationScorerFactory : public ScorerFactory {
 public:
  /**
   * Default constructor
   */
  CrossValidationScorerFactory() = default;

  /**
   * Create an instance of a #sgpp::datadriven::CrossValidation object based on the configuration
   * @param parser Instance of #sgpp::datadriven::DataMiningConfigParser that reads the required
   * data from the config file.
   * @return Fully configured instance of a  #sgpp::datadriven::CrossValidation object.
   */
  Scorer *buildScorer(const DataMiningConfigParser &parser) const override;
};

} /* namespace datadriven */
} /* namespace sgpp */
