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
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationBit.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/parameters/ContinuousParameter.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/parameters/DiscreteParameter.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Concrete factory to build instances of #sgpp::datadriven::ModelFittingDensityEstimation
 */
class DensityEstimationFitterFactory : public FitterFactory {
 public:
  /**
   * Default constructor
   */
  DensityEstimationFitterFactory(DataMiningConfigParser &parser);

  /**
   * Assemble a #sgpp::datadriven::ModelFittingDensityEstimation object based on the configuration
   * determined by a previous set_() call.
   * @return Fully configured instance of a  #sgpp::datadriven::ModelFittingDensityEstimation object.
   */
  ModelFittingBase *buildFitter() override;

  /**
   * Outputs information about the current hyperparameter configuration.
   * @return String to print to console or file containing values of manipulated hyperparameters.
   */
  std::string printConfig() override;

 protected:
  /**
   * Configuration for all parameters that are not optimized
   */
  FitterConfigurationDensityEstimation baseConfig;
  /**
   * Container for GridTypes specifically for the basis function hyperparameter
   */
  std::vector<base::GridType> basisFunctions;

};

} /* namespace datadriven */
} /* namespace sgpp */
