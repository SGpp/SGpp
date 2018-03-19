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
  DiscreteParameter() = default;

	DiscreteParameter(std::string&& name, int min, int max);

  DiscreteParameter(int nBits, std::string &name, int min, int max)
          :HyperParameter(nBits, name), min(min), max(max){}
	//~DiscreteParameter();

  int getValue();
	int getNOptions();
	void setBO(int option);
	void setHarmonica() override;

protected:
	int min;
	int max;
	int value = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif /* DISCRETEPARAMETER_HPP_ */
