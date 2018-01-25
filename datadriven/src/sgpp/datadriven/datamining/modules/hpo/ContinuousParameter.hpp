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
	ContinuousParameter(double min, double max):min(min),max(max){}
	virtual ~ContinuousParameter();
	double getValue(int* configID);


protected:
	double min;
	double max;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif /* CONTINUOUSPARAMETER_HPP_ */
