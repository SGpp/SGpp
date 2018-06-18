/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DiscreteParameter.hpp
 *
 *  Created on: Jan 25, 2018
 *      Author: Eric Koepke
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
   * @param option value of the parameter from 0
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
