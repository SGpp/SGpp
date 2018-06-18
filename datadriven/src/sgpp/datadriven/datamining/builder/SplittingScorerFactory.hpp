/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SplittingScorerFactory.hpp
 *
 * Created on: Oct 11, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/builder/ScorerFactory.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Concrete factory to build an instance of #sgpp::datadriven::SplittingScorer
 */
class SplittingScorerFactory : public ScorerFactory {
 public:
  /**
   * Default constructor
   */
  SplittingScorerFactory() = default;

  /**
   * Create an instance of a #sgpp::datadriven::SplittingScorer object based on the configuration
   * @param parser Instance of #sgpp::datadriven::DataMiningConfigParser that reads the required
   * data from the config file.
   * @return Fully configured instance of a  #sgpp::datadriven::SplittingScorer object.
   */
  Scorer *buildScorer(const DataMiningConfigParser &parser) const override;
};
} /* namespace datadriven */
} /* namespace sgpp */
