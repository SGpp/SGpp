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

#include <string>

#include <sgpp/datadriven/datamining/modules/hpo/parameters/HyperParameter.hpp>

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
   * @param minv minimum value of the hyperparameter during optimization
   * @param maxv maximum value of the hyperparameter during optimization
   */
  DiscreteParameter(std::string && name, int minv, int maxv);

  /**
   * Constructor with custom number of bits
   * @param nBits number of bits for representation in harmonica
   * @param name name of the hyperparameter
   * @param minv minimum value of the hyperparameter during optimization
   * @param maxv maximum value of the hyperparameter during optimization
   */
  DiscreteParameter(int nBits, std::string && name, int minv, int maxv)
      : HyperParameter(nBits, name), minv(minv), maxv(maxv) {}

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
  int minv;
  /**
   * maximum value of the hyperparameter during optimization
   */
  int maxv;
  /**
   * current value of the hyperparameter
   */
  int value = 0;
};
} /* namespace datadriven */
} /* namespace sgpp */
#endif /* DISCRETEPARAMETER_HPP_ */
