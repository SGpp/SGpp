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
#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationBit.hpp>


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
  LeastSquaresRegressionFitterFactory(DataMiningConfigParser parser);

  /**
   * Assemble a #sgpp::datadriven::ModelFittingBase object based on the configuration
   * @param parser Instance of #sgpp::datadriven::DataMiningConfigParser that reads the required
   * data from the config file.
   * @return Fully configured instance of a  #sgpp::datadriven::ModelFittingBase object.
   */
  ModelFittingBase* buildFitter(int configID, int row, DataMatrix paritymatrix) override;
  
  FitterConfiguration* buildConfig() const override;
  
  int buildParity() override;

  void addConstraint(int idx, int bias) override;

 protected:
  std::list<ConfigurationBit> configBits;
  FitterConfigurationLeastSquares baseConfig;
  std::map<std::string,ContinuousParameter> conpar;
  std::map<std::string,DiscreteParameter> dispar;
  std::vector<std::list<ConfigurationBit*> > parityrow;
};

} /* namespace datadriven */
} /* namespace sgpp */
