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

#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>


namespace sgpp {
namespace datadriven {




void FitterFactory::setHarmonica()  {

  for(auto& pair: conpar){
    pair.second.setHarmonica();
  }
  for(auto& pair: dispar){
    pair.second.setHarmonica();
  }
  for(auto& pair: catpar){
    pair.second.setHarmonica();
  }
  // std::cout<<"Run mark 2.1"<<std::endl;
}




void FitterFactory::setBO(BOConfig* config){ //EDIT: old one destroyed
  int i = 0;
  for(auto& pair: dispar){
    pair.second.setBO(config->getDisc(i));
    i++;
  }
  i = 0;
  for(auto& pair: catpar){
    pair.second.setBO(config->getCat(i));
    i++;
  }
  i = 0;
  for(auto& pair: conpar){
    pair.second.setBO(config->getCont(i));
    i++;
  }
}

void FitterFactory::getConfigBits(std::vector<ConfigurationBit*>& configBits) {
  for(auto& pair: conpar){
    pair.second.makeConfigBits(configBits);
  }
  for(auto& pair: dispar){
    pair.second.makeConfigBits(configBits);
  }
  for(auto& pair: catpar){
    pair.second.makeConfigBits(configBits);
  }

}

BOConfig FitterFactory::getBOConfig() {
  for(auto& pair: dispar){
    discOptions.push_back(pair.second.getNOptions());
  }
  for(auto& pair: catpar){
    catOptions.push_back(pair.second.getNOptions());
  }
  return BOConfig(&discOptions, &catOptions, conpar.size());
}

} /* namespace datadriven */
} /* namespace sgpp */
