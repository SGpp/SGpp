/*
 * ContinuousParameter.hpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#ifndef CONTINUOUSPARAMETER_HPP_
#define CONTINUOUSPARAMETER_HPP_

#include "HyperParameter.hpp"

namespace sgpp {
namespace datadriven {

class ContinuousParameter: public sgpp::datadriven::HyperParameter {
public:
  ContinuousParameter() = default;
	ContinuousParameter(int nBits, std::string&& name, double min, double max)
          :HyperParameter(nBits, name), min(min), max(max){}
	// ~ContinuousParameter();
	virtual double getValue();
	void setHarmonica() override;
	void setBO(double interval);


protected:
	double min;
	double max;
	double value = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif /* CONTINUOUSPARAMETER_HPP_ */
