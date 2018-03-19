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


int ConfigurationBit::evaluate(){
  if(bVisited){
	  //std::cout<<"Bit visited!"<<std::endl;
    return value;
  }
  bVisited = true;
  //for each constraint
  for(auto &constraint : constraints){
    if(constraint->getConfigBits().size() == 1){
      value = constraint->getBias();
      return value;
    }
  }
  for(auto &constraint : constraints){
    int tmp = constraint->getBias();
    for(auto &bit : constraint->getConfigBits()){
    	if(this!=bit){
    		tmp = tmp * bit->evaluate();
    	}
    }
    if(tmp != 0){
      value = tmp;
      return value;
    }
  }
  // pull free bit
  std::cout<<"Error: Bit evaluation failed."<<std::endl;
  std::cout<<name<<std::endl;
  for(auto &constraint : constraints){
    std::cout<<"Bias: "<<constraint->getBias()<<std::endl;
    for(auto &bit : constraint->getConfigBits()){
      std::cout<<"Bit "<<bit->name<<": "<<bit->evaluate()<<std::endl;
    }
  }
  // std::cout<<"input after:"<<*input<<std::endl;
  return value;
}

int ConfigurationBit::fixFreeBits(std::vector<ConfigurationBit*> &freeBits){
  if(bVisited){
    return value;
  }
  bVisited = true;
  //for each constraint
  for(auto &constraint : constraints){
    if(constraint->getConfigBits().size() == 1){
      value = constraint->getBias();
      // std::cout<<"Constraint Single Bias: "<<value<<std::endl;
      return value;
    }
  }
  for(auto &constraint : constraints){
    int tmp = constraint->getBias();
    // std::cout<<"Constraint Bias: "<<tmp<<std::endl;
    for(auto &bit : constraint->getConfigBits()){
    	if(this!=bit){
    		tmp = tmp * bit->fixFreeBits(freeBits);
    		// std::cout<<"temp: "<<tmp<<std::endl;
    	}
    }
    if(tmp != 0){
        // std::cout<<"Constraint Resolved!"<<std::endl;
      value = tmp;
      return value;
    }
  }
  // pull free bit
  value = 1;
  freeBits.push_back(this);
  return value;
}

void ConfigurationBit::addConstraint(ConfigurationRestriction* constraint){
  constraints.push_back(constraint);
}

void ConfigurationBit::removeLastConstraint(){
  constraints.pop_back();
}

void ConfigurationBit::reset(){
  value = 0;
  bVisited = false;
}

void ConfigurationBit::setBit(int* input){
  value = (*input&1)*2-1;
  *input = (*input>>1);
  bVisited = true;
}

bool ConfigurationBit::checkConstraints(){
	 for(auto &constraint : constraints){
	    int tmp = constraint->getBias();
	    // std::cout<<"Constraint Bias: "<<tmp<<std::endl;
	    for(auto &bit : constraint->getConfigBits()){
	    	if(this!=bit){
	    		tmp = tmp * bit->evaluate();
	    		// std::cout<<"temp: "<<tmp<<std::endl;
	    	}
	    }
	    if(tmp == 0){
	        std::cout<<"Error: Unset Bit on Constraint Evaluation"<<std::endl;
	    }
	    if(value != tmp){
	    	return false;
	    }

	  }
return true;
}

bool ConfigurationBit::setValue(int input) {
  for(auto &constraint : constraints) {
    constraint->reduceOpenBits();
  }
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
  return true;
}

int ConfigurationBit::getValue() {
  return value;
}


}  // namespace datadriven
}  // namespace sgpp
