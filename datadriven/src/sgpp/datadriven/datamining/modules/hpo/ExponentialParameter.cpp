/*
 * CategoricalParameter.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#include "ExponentialParameter.hpp"

#include <math.h>
#include <cmath>

namespace sgpp {
namespace datadriven {




double ExponentialParameter::getValue(){
	return pow(10,value);
}


} /* namespace datadriven */
} /* namespace sgpp */
