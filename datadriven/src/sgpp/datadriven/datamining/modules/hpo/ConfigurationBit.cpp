// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationBit.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>
#include <list>
#include <iostream>

namespace sgpp {
namespace datadriven {



void ConfigurationBit::addConstraint(ConfigurationRestriction* constraint){
  constraints.push_back(constraint);
}

void ConfigurationBit::removeLastConstraint(){
  constraints.pop_back();
}

void ConfigurationBit::reset(){
  value = 0;
}



void ConfigurationBit::setValue(int input) {
  value = input;


  for(auto &constraint : constraints) {
   // std::cout<<"Test Point 178!"<<std::endl; //EDIT: throw exception
    constraint->reduceOpenBits();
  }
 // std::cout<<"Test Point 188!"<<std::endl; //EDIT: throw exception

  /*for(auto &constraint : constraints){
    int nOpen = constraint->getOpenBits();
    if(nOpen == 0){
      if(!constraint->check()){
        return false;
      }
    }else if(nOpen == 1){
      if(!constraint->resolve()){
        return false;
      }
    }
  }*/
}

int ConfigurationBit::getValue() {
  return value;
}

std::string ConfigurationBit::getName() {
  return name;
}


}  // namespace datadriven
}  // namespace sgpp
