/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
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
   * @param logscale whether this parameter operates on a logscale
   */
  ContinuousParameter(size_t nBits,
                      std::string &&name,
                      double min,
                      double max,
                      bool logscale = false)
      : HyperParameter(nBits, name), min(min), max(max), logscale(logscale) {}
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
  bool logscale = false;
};
} /* namespace datadriven */
} /* namespace sgpp */
#endif /* CONTINUOUSPARAMETER_HPP_ */
