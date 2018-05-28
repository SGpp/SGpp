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

/**
 * concrete class for hyperparameter with discrete values
 */
class DiscreteParameter : public sgpp::datadriven::HyperParameter {
 public:
  /**
   * Default Constructor
   */
  DiscreteParameter() = default;

  /**
   * Normal constructor, number of bits calculated automatically
   * @param name name of the hyperparameter
   * @param min minimum value of the hyperparameter during optimization
   * @param max maximum value of the hyperparameter during optimization
   */
  DiscreteParameter(std::string &&name, int min, int max);

  /**
   * Constructor with custom number of bits
   * @param nBits number of bits for representation in harmonica
   * @param name name of the hyperparameter
   * @param min minimum value of the hyperparameter during optimization
   * @param max maximum value of the hyperparameter during optimization
   */
  DiscreteParameter(int nBits, std::string &name, int min, int max)
      : HyperParameter(nBits, name), min(min), max(max) {}

  /**
  * Retrieve the current value of the hyperparameter
  * @return value of the hyperparameter
  */
  int getValue();

  /**
   * Retrieve the number of options for this parameter
   * @return
   */
  int getNOptions();
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
  void setBO(int option);

 protected:
  /**
   * minimum value of the hyperparameter during optimization
   */
  int min;
  /**
   * maximum value of the hyperparameter during optimization
   */
  int max;
  /**
   * current value of the hyperparameter
   */
  int value = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */
#endif /* DISCRETEPARAMETER_HPP_ */
