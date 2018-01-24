// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationBit.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {


int evaluate(int* input){
  if(bVisited){
    return value;
  }
  bVisited = true;
  //for each constraint
  for(const auto &constraint : constraints){
    if(constraint.getConfigBits().size() == 0){
      value = constraint.getBias();
      return value;
    }
  }
  for(const auto &constraint : constraints){
    int tmp = constraint.getBias();
    for(const auto &bit : constraint.getConfigBits()){
      tmp = tmp * bit.evaluate();
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


void addConstraint(ConfigurationRestriction* constraint){
  constraints.append(constraint);
}

void reset(){
  value = 0;
  bVisited = false;
}


}  // namespace datadriven
}  // namespace sgpp
