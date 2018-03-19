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

void FitterFactory::getBOspace(int* nCont, std::vector<int>& nOptions){

  *nCont = conpar.size();
  for(auto& pair: dispar){
    nOptions.push_back(pair.second.getNOptions());
  }
  for(auto& pair: catpar){
    nOptions.push_back(pair.second.getNOptions());
  }

} //EDIT: add categorical parameters

void FitterFactory::setBO(base::DataVector& cont, std::vector<int>& disc){
  int i = 0;
  for(auto& pair: dispar){
    pair.second.setBO(disc[i]);
    i++;
  }
  for(auto& pair: catpar){
    pair.second.setBO(disc[i]);
    i++;
  }
  i = 0;
  for(auto& pair: conpar){
    pair.second.setBO(cont[i]);
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

} /* namespace datadriven */
} /* namespace sgpp */
