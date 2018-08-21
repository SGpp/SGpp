/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterConfigurationLeastSquares.hpp
 *
 * Created on: Oct 10, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Configuration for fitter scenarios using least squares optimization.
 */
class FitterConfigurationLeastSquares : public FitterConfiguration {
 public:
  FitterConfigurationLeastSquares() = default;

  FitterConfiguration *clone() const override;

  void setupDefaults() override;

  /**
   * First setup default values, then read new input values from configuration file.
   * @param parser the parsed configuration file.
   */
  void readParams(const DataMiningConfigParser &parser) override;
};
} /* namespace datadriven */
} /* namespace sgpp */
