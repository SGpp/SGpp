/*
 * ContinuousParameter.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#include "ContinuousParameter.hpp"

namespace sgpp {
namespace datadriven {



void ContinuousParameter::setHarmonica() {
  double v = 0;
  double m = 1;
  for(auto &bit : bits){
    v = v + m* bit->evaluate();
    m = m * 2;
  }
  value = min+((max-min)*(1+v/(m-1.0))/2);
}

void ContinuousParameter::setBO(double interval) {
  value = (max-min)*interval + min;
}

double ContinuousParameter::getValue(){
  return value;
}

} /* namespace datadriven */
} /* namespace sgpp */
