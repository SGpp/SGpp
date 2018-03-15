/*
 * CategoricalParameter.hpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#ifndef CategoricalParameter_HPP_
#define CategoricalParameter_HPP_

#include "HyperParameter.hpp"

namespace sgpp {
namespace datadriven {

template <class T>
class CategoricalParameter: public sgpp::datadriven::HyperParameter {
public:
	CategoricalParameter(int min, int max):min(min),max(max){}
	//~CategoricalParameter();
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
#endif /* CategoricalParameter_HPP_ */
