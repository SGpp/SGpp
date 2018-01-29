/*
 * ContinuousParameter.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#include "ContinuousParameter.hpp"

namespace sgpp {
namespace datadriven {





double ContinuousParameter::getValue(int* configID){
	double v = 0;
	double m = 1;
  for(auto &bit : bits){
    v = v + m* bit->evaluate(configID);
    m = m * 2;
  }
  return min+((max-min)*(1+v/(m-1.0))/2);
}

} /* namespace datadriven */
} /* namespace sgpp */
