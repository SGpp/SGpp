/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CombiConfigurator.hpp
 *
 *  Created on: Jan 29, 2019
 *      Author: nico
 */
#pragma once
#ifdef USE_PYTHON_EMBEDDING
#include <python3.6/Python.h>
#endif
#include <iostream>
#include <vector>

namespace sgpp {
namespace datadriven {

using std::vector;

struct combiConfig {
  double coef;
  vector<size_t> levels;
};

class CombiConfigurator {
 public:
  CombiConfigurator();

  /**
   * Initializes an adaptive combination scheme in SpACE
   * @param dim Dimension of the sparse grid
   * @param level of the initial sparse grid
   */
  void initAdaptiveScheme(size_t dim, size_t level);

  /**
   * Delivers the current set of component configs
   * @param vec reference for storing the vector of component configs
   */
  void getCombiScheme(vector<combiConfig> &vec);

  /**
   * Checks if a component is refineable
   * @param levelvec config of the component
   * @return returns true when the component is refinable
   */
  bool isRefinable(combiConfig levelvec);

  /**
   * Refines a single component and therefore changes the set of components
   * @param levelvec config of the component
   */
  void refineComponent(combiConfig levelvec);

 private:
  /**
   * Initializes the Python Interpreter
   */
  void initializePython();

  /**
   * Finalizes the Python Interpreter
   */
  void finalizePython();

#ifdef USE_PYTHON_EMBEDDING
  /**
   * Converts a components config in SpACE to one in SG++
   * @param *pValue PyObject that corrensponds to the config in SpACE
   * @return returns a combiConfig
   */
  inline combiConfig combiConfFromPyObj(PyObject *pValue);

  /**
   * Convers a components config in SG++ to one in SpACE
   * @param pair combiConfig from SG++
   * @return returns a components config in SpACE as PyObject
   */
  inline PyObject *combiConfAsPyObj(combiConfig pair);

  /// Instance of the SpACE CombiScheme
  PyObject *combischeme;
#endif
};
} /* namespace datadriven */
} /* namespace sgpp */
