/*
 * DiscreteParameter.hpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#ifndef DISCRETEPARAMETER_HPP_
#define DISCRETEPARAMETER_HPP_

#include "HyperParameter.hpp"

namespace sgpp {
namespace datadriven {

class DiscreteParameter: public sgpp::datadriven::HyperParameter {
public:
	DiscreteParameter(int min, int max):min(min),max(max){}
	virtual ~DiscreteParameter();
	std::list<ConfigurationBit> makeConfigBits();
	int getValue(int* configID);

protected:
	int min;
	int max;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif /* DISCRETEPARAMETER_HPP_ */
