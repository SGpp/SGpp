/*
 * HyperParameter.cpp
 *
 *  Created on: 24.01.2018
 *      Author: Eric
 */

#include <iostream>
#include "HyperParameter.hpp"

namespace sgpp {
namespace datadriven {



void HyperParameter::makeConfigBits(std::vector<ConfigurationBit*>& configBits){
  bits.reserve(nBits);
	for(int i=0;i<nBits;i++){
    std::string bitName = name + std::to_string(i);
    bits.emplace_back(bitName);
    configBits.push_back(&bits[i]);
  }

}

} /* namespace datadriven */
} /* namespace sgpp */
