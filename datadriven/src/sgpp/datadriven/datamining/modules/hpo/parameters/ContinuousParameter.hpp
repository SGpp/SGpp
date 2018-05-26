/*
 * ContinuousParameter.hpp
 *
 *  Created on: Jan 25, 2018
 *      Author: Eric Koepke
 */

#ifndef CONTINUOUSPARAMETER_HPP_
#define CONTINUOUSPARAMETER_HPP_

#include "HyperParameter.hpp"

namespace sgpp {
namespace datadriven {

/**
 * Concrete class for hyperparameter with continuous values
 */
class ContinuousParameter: public sgpp::datadriven::HyperParameter {
public:
  /**
   * Default constructor
   */
  ContinuousParameter() = default;
  /**
   * Normal Constructor
   * @param nBits 
   * @param name
   * @param min
   * @param max
   */
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
