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
	//~DiscreteParameter();
	void makeConfigBits(std::list<std::unique_ptr<ConfigurationBit>>& allbits);
	int getValue();
	int getNOptions();
	void setBO(int option);
	void setHarmonica() override;

protected:
	int min;
	int max;
	int value;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif /* DISCRETEPARAMETER_HPP_ */
