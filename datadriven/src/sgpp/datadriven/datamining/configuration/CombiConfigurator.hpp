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
  bool refineBlock(vector<combiConfig> &vec, size_t index);

  void test(vector<combiConfig> &vec);

 private:
  void Initialize();
  void Finalize();

  inline static PyObject *pName, *pModule, *pFunc;
  inline static bool initialized = false;
};
} /* namespace datadriven */
} /* namespace sgpp */
