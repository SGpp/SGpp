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
	DiscreteParameter(std::string& name, int min, int max);
  DiscreteParameter(std::string &name, int min, int max, int nBits)
          : HyperParameter(name), min(min), max(max), nBits(nBits) {}
	//~DiscreteParameter();
	void makeConfigBits(std::list<std::unique_ptr<ConfigurationBit>>& allbits);

  virtual int getValue();
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
