/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * FitterFactory.hpp
 *
 * Created on: Dec 17, 2017
 *     Author: Eric Koepke
 */

#pragma once

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>


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
   * @param parser Instance of #sgpp::datadriven::DataMiningConfigParser that reads the required
   * data from the config file.
   * @return Fully configured instance of a  #sgpp::datadriven::ModelFittingBase object.
   */
  
  virtual FitterConfiguration* buildConfig() const = 0;

  virtual ModelFittingBase* buildFitter(int configID, int row, DataMatrix &paritymatrix);

  virtual int buildParity();

  virtual int addConstraint(int idx, int bias);

};
} /* namespace datadriven */
} /* namespace sgpp */
