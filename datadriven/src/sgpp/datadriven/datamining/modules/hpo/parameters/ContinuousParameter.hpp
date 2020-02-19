// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CONTINUOUSPARAMETER_HPP_
#define CONTINUOUSPARAMETER_HPP_

#include <sgpp/datadriven/datamining/modules/hpo/parameters/HyperParameter.hpp>

#include <string>

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
   * @param minv minimum value of the hyperparameter during optimization
   * @param maxv maximum value of the hyperparameter during optimization
   * @param logscale whether this parameter operates on a logscale
   */
  ContinuousParameter(size_t nBits,
                      std::string && name,
                      double minv,
                      double maxv,
                      bool logscale = false)
      : HyperParameter(nBits, name), minv(minv), maxv(maxv), logscale(logscale) {}
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
  double minv;
  /**
   * maximum value of the hyperparameter during optimization
   */
  double maxv;
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
