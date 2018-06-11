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
class ContinuousParameter : public sgpp::datadriven::HyperParameter {
 public:
  /**
   * Default constructor
   */
  ContinuousParameter() = default;
  /**
   * Normal Constructor
   * @param nBits number of bits for representation in harmonica
   * @param name name of the hyperparameter
   * @param min minimum value of the hyperparameter during optimization
   * @param max maximum value of the hyperparameter during optimization
   */
  ContinuousParameter(int nBits, std::string &&name, double min, double max)
      : HyperParameter(nBits, name), min(min), max(max) {}
  // ~ContinuousParameter();
  /**
   * Retrieve the current value of the hyperparameter
   * @return value of the hyperparameter
   */
  virtual double getValue();
  /**
   * adjust the current value of the hyperparameter according to the bit
   * configuration by harmonica
   */
  void setHarmonica() override;
  /**
   * adjust the current value of the hyperparameter according to the (normalized)
   * input
   * @param interval (normalized) value of the hyperparameter
   */
  void setBO(double interval);

 protected:
  /**
   * minimum value of the hyperparameter during optimization
   */
  double min;
  /**
   * maximum value of the hyperparameter during optimization
   */
  double max;
  /**
   * current value of the hyperparameter
   */
  double value = 0;

  /**
   * whether the parameter is on a log-scale
   */
  bool logscale = False;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif /* CONTINUOUSPARAMETER_HPP_ */
