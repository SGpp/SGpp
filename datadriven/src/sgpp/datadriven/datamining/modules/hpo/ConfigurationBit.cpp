// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationBit.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>
#include <list>

namespace sgpp {
namespace datadriven {


int ConfigurationBit::evaluate(int* input){
  if(bVisited){
    return value;
  }
  bVisited = true;
  //for each constraint
  for(auto &constraint : constraints){
    if(constraint.getConfigBits().size() == 0){
      value = constraint.getBias();
      return value;
    }
  }
  for(auto &constraint : constraints){
    int tmp = constraint.getBias();
    for(auto &bit : constraint.getConfigBits()){
      tmp = tmp * bit.evaluate(input);
    }
    if(tmp != 0){
      value = tmp;
      return value;
    }
  }
  // pull free bit
  value = (*input&1)*2-1;
  *input = (input>>1);
  return value;
}


void ConfigurationBit::addConstraint(ConfigurationRestriction* constraint){
  constraints.push_back(*constraint);
}

void ConfigurationBit::reset(){
  value = 0;
  bVisited = false;
}


}  // namespace datadriven
}  // namespace sgpp
