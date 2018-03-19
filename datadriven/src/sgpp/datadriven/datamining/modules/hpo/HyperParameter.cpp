/*
 * HyperParameter.cpp
 *
 *  Created on: 24.01.2018
 *      Author: Eric
 */

#include "HyperParameter.hpp"

namespace sgpp {
namespace datadriven {

HyperParameter::HyperParameter() {
	// TODO Auto-generated constructor stub

}

HyperParameter::~HyperParameter() {
	// TODO Auto-generated destructor stub
}

void HyperParameter::makeConfigBits(std::vector<ConfigurationBit*>& configBits){
	for(int i=0;i<nBits;i++){
    bits.push_back(std::make_unique<ConfigurationBit>(name + std::to_string(i)));
    configBits.push_back(bits[i].get());
	}
}

} /* namespace datadriven */
} /* namespace sgpp */
