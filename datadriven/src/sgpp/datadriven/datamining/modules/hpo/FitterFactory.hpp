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
#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationBit.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ContinuousParameter.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/DiscreteParameter.hpp>


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

  virtual ModelFittingBase* buildFitter() = 0;

  int buildParity();

  int addConstraint(int idx, int bias);

  virtual void printConfig() = 0;

  void setHarmonica(int configID, int row, DataMatrix &paritymatrix);

  void getBOspace(int* nCont, std::vector<int>& nOptions); //EDIT: add categorical parameters

  void setBO(base::DataVector& cont, std::vector<int>& disc);

protected:
  std::list<std::unique_ptr<ConfigurationBit>> configBits;
  std::map<std::string,ContinuousParameter*> conpar;
  std::map<std::string,DiscreteParameter*> dispar;
  std::vector<std::list<ConfigurationBit*> > parityrow;
  std::vector<ConfigurationBit*> freeBits;


  };
} /* namespace datadriven */
} /* namespace sgpp */
