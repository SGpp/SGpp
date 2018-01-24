/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * LeastSquaresRegressionFitterFactory.cpp
 *
 *  Created on:	17.12.2017
 *      Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/builder/LeastSquaresRegressionFitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquares.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

namespace sgpp {
namespace datadriven {

LeastSquaresRegressionFitterFactory::LeastSquaresRegressionFitterFactory()
  :configBits(){
  // build ConfigurationBits (constructor)
  for(int i=0;i<12;i++){
    configBits.append(new ConfigurationBit());
  }
}

ModelFittingBase* LeastSquaresRegressionFitterFactory::buildFitter(int configID) const {
  // fix ConfigurationBits according to constraints
  for(auto &bit : configBits){
    bit.reset();
  }
  for(auto &bit : configBits){
    bit.evaluate(&configID);
  }

  // build config
  
  // return model and ConfigurationBits in vector/matrix
  
  return new ModelFittingLeastSquares(*static_cast<FitterConfigurationLeastSquares*>(config));
}
FitterConfiguration* LeastSquaresRegressionFitterFactory::buildConfig() const {
  FitterConfigurationLeastSquares* config = new FitterConfigurationLeastSquares();
  DataMiningConfigParser parser("TrainWithTestingExample.json");
  config->readParams(parser);
  return config;
}
} /* namespace datadriven */
} /* namespace sgpp */
