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

std::list<ConfigurationBit> HyperParameter::makeConfigBits(int nBits){
	for(int i=0;i<nBits;i++){
		bits.push_back(*(new ConfigurationBit()));
	}
	return bits;
}

} /* namespace datadriven */
} /* namespace sgpp */
