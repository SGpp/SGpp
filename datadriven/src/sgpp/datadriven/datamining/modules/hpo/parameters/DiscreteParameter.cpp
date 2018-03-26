/*
 * DiscreteParameter.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#include "DiscreteParameter.hpp"

#include <math.h>
#include <cmath>

namespace sgpp {
namespace datadriven {


DiscreteParameter::DiscreteParameter(std::string&& name, int min, int max)
        :HyperParameter(0, name), min(min), max(max){
  int i = 1;
  int c = 2;
  while(c < max-min+1){
    i++;
    c=c*2;
  }
  nBits = i;
}



int DiscreteParameter::getValue(){
	return value;
}

int DiscreteParameter::getNOptions() {
  return max-min+1;
}

void DiscreteParameter::setBO(int option){
  value = min+option;
}
void DiscreteParameter::setHarmonica(){
  //if(std::pow(2,bits.size()) <= max-min+1){
    double v = 0;
    double m = 1;
    for(auto bit : bits){
      v = v + m* bit.getValue();
      m = m * 2;
    }
    value = lround(min+((max-min)*(1.0+v/(m-1.0))/2.0));
  //}else{ // if(std::pow(2,bits.size()) > max-min+1)
    // EDIT: try different ways here later, like above and with max lenght c
    /* int c = std::pow(2,bits.size()-1);
    int k = 0;
    while(k+c/2+max-min+1 <= std::pow(2,bits.size())){
      c = c/2;
      k = k + c;
    }
    int c = 1;
    while(c < std::pow(2,bits.size())-(max-min+1)){
      c = c*2;
    }
    // k = k - c;
    int v = min;
    int m = 1;
    // make sure this doesn't break it
    ConfigurationBit* last = bits.back();
    bits.pop_back();
    for(auto bit : bits){
      v = v + m*(bit->evaluate()+1)/2;
      m = m * 2;
    }
    bits.push_back(last);
    int nv = v + c*(last->evaluate()+1)/2;
    if(nv >= min+std::pow(2,bits.size()-1) && nv <= max){
      v = nv;
    }
    value = v;
  }*/
}

} /* namespace datadriven */
} /* namespace sgpp */
