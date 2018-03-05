/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterConfigurationDensityEstimation.hpp
 *
 * Created on: Jan 02, 2018
 *     Author: Kilian Röhner
 */

#pragma once

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Configuration for fitter scenarios using density estimation.
 */
class FitterConfigurationDensityEstimation : public FitterConfiguration {
 public:
  FitterConfigurationDensityEstimation() = default;

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
