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

void HyperParameter::makeConfigBits(int nBits, std::list<std::unique_ptr<ConfigurationBit>>& allbits){
	for(int i=0;i<nBits;i++){
		ConfigurationBit* n = new ConfigurationBit();
		bits.push_back(n);
		allbits.push_back(std::unique_ptr<ConfigurationBit>(n));
	}
}

} /* namespace datadriven */
} /* namespace sgpp */
