/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * LeastSquaresRegressionFitterFactory.hpp
 *
 *  Created on:	17.12.2017
 *      Author: Eric Koepke
 */

#pragma once

#include <sgpp/datadriven/datamining/builder/FitterFactory.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Concrete factory to build an instance of #sgpp::datadriven::ModelFittingBase
 */
class LeastSquaresRegressionFitterFactory : public FitterFactory {
 public:
  /**
   * Default constructor
   */
  LeastSquaresRegressionFitterFactory();

  /**
   * Assemble a #sgpp::datadriven::ModelFittingBase object based on the configuration
   * @param parser Instance of #sgpp::datadriven::DataMiningConfigParser that reads the required
   * data from the config file.
   * @return Fully configured instance of a  #sgpp::datadriven::ModelFittingBase object.
   */
  ModelFittingBase* buildFitter(FitterConfiguration* config) const override;
  
  FitterConfiguration* buildConfig() const override;
  
 protected:
  std::list<ConfigurationBit> configBits;
};

} /* namespace datadriven */
} /* namespace sgpp */
