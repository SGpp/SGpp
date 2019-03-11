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

#include <python3.7/Python.h>
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

  void initAdaptiveScheme(size_t dim, size_t level);

  void getCombiScheme(vector<combiConfig> &vec);

  bool isRefinable(combiConfig levelvec);

  void refineComponent(combiConfig levelvec);

  /**
   * This returns the set of level vectors with the corresponding coefficients
   * for the standard combination technique
   * @param vec reference where to store the combiConfigSet
   * @oaram d dimension of the sparse grid
   * @param l level of the sparse grid
   * @return vector of combiConfig
   */
  void getStandardCombi(vector<combiConfig> &vec, size_t d, size_t l);

  /**
   * @param vec reference where to store the combiConfigSet
   * @param index to the block that needs to be refined
   * @return new vector of combiConfig
   * @return true if the block could be refined
   */

  void test(vector<combiConfig> &vec);

 private:
  void initializePython();
  void finalizePython();

  inline combiConfig combiConfFromPyObj(PyObject *pValue);
  inline PyObject *combiConfAsPyObj(combiConfig pair);

  PyObject *combischeme;

  inline static PyObject *pModule;
  inline static bool initialized = false;
};
} /* namespace datadriven */
} /* namespace sgpp */
