/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterConfigurationClassification.hpp
 *
 *  Created on: Jul 2, 2018
 *      Author: dominik
 */

#pragma once

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Configuration for fitter scenarios using classification
 */
class FitterConfigurationClassification : public FitterConfigurationDensityEstimation {
 public:
  FitterConfigurationClassification() = default;

  FitterConfiguration* clone() const override;

  void setupDefaults() override;

  /**
   * First setup default values, then read new input values from configuration file.
   * @param parser the parsed configuration file.
   */
  void readParams(const DataMiningConfigParser& parser) override;
};

} /* namespace datadriven */
} /* namespace sgpp */

