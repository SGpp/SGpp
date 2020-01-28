// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationBit.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/parameters/ContinuousParameter.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/parameters/DiscreteParameter.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BOConfig.hpp>

#include <vector>
#include <map>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Abstract factory to build all kinds of fitters/models based on a given configuration.
 */
class FitterFactory {
 public:
  /**
   * Default constructor
   */
  FitterFactory() = default;

  /**
   * Virtual destructor
   */
  virtual ~FitterFactory() = default;

  /**
   * Assemble a #sgpp::datadriven::ModelFittingBase object based on the configuration
   * determined by a previous set_() call.
   * @return Fully configured instance of a  #sgpp::datadriven::ModelFittingBase object.
   */
  virtual ModelFittingBase *buildFitter() = 0;

  /**
   * Outputs information about the current hyperparameter configuration.
   * @return String to print to console or file containing values of manipulated hyperparameters.
   */
  std::string printConfig();

  /**
   * Outputs the names of the hyperparameters
   * @return string to pringt ot console or file containing comma seperated names
   */
  std::string printHeadline();

  /**
   * Setup connection to hyperparameter classes for modifying them through boolean represenations
   * @param configBits reference to the boolean "configBits"
   */
  virtual void getConfigBits(std::vector<ConfigurationBit *> &configBits);

  /**
   * Adjusts the current hyperparameter configuration according to the configBits
   */
  void setHarmonica();

  /**
   * Gets a compact representation of the hyperparameter configuration space.
   * @return Object of type BOConfig that can be cloned to create different
   * hyperparameter configurations
   */
  BOConfig getBOConfig();

  /**
   * Adjusts current hyperparameter configuration according to the input
   * @param config compact representation of hyperparameters
   */
  void setBO(BOConfig &config);

 protected:
  /**
   * map to store hyperparameters defined on a continuous domain
   */
  std::map<std::string, ContinuousParameter> conpar;
  /**
   * map to store hyperparameters defined on a discrete domain
   */
  std::map<std::string, DiscreteParameter> dispar;
  /**
    * map to store hyperparameters defined on a discrete domain
    * without inherent ordering (categorical)
    */
  std::map<std::string, DiscreteParameter> catpar;
  /**
   * number of options for all discrete parameters
   */
  std::vector<int> discOptions;
  /**
   * numer of options for all categorical parameters
   */
  std::vector<int> catOptions;

  /**
   * Container for GridTypes specifically for the basis function hyperparameter
   */
  std::vector<base::GridType> basisFunctions;
};
} /* namespace datadriven */
} /* namespace sgpp */
