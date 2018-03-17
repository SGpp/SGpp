/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DensityEstimationFitterFactory.hpp
 *
 *  Created on:	17.12.2017
 *      Author: Eric Koepke
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationBit.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ContinuousParameter.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/DiscreteParameter.hpp>




namespace sgpp {
namespace datadriven {

/**
 * Concrete factory to build an instance of #sgpp::datadriven::ModelFittingBase
 */
class DensityEstimationFitterFactory : public FitterFactory {
 public:
  /**
   * Default constructor
   */
  DensityEstimationFitterFactory(DataMiningConfigParser& parser);

  /**
   * Assemble a #sgpp::datadriven::ModelFittingBase object based on the configuration
   * @param parser Instance of #sgpp::datadriven::DataMiningConfigParser that reads the required
   * data from the config file.
   * @return Fully configured instance of a  #sgpp::datadriven::ModelFittingBase object.
   */
  ModelFittingBase* buildFitter() override;

  void printConfig() override;


protected:
  FitterConfigurationLeastSquares baseConfig;

};

} /* namespace datadriven */
} /* namespace sgpp */
