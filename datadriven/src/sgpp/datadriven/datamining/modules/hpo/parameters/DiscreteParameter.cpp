/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DiscreteParameter.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: Eric Koepke
 */

#include <sgpp/datadriven/datamining/modules/hpo/parameters/DiscreteParameter.hpp>

#include <cmath>
#include <string>

namespace sgpp {
namespace datadriven {

DiscreteParameter::DiscreteParameter(std::string && name, int minv, int maxv)
    : HyperParameter(0, name), minv(minv), maxv(maxv) {
  size_t i = 1;
  int c = 2;
  while (c < maxv - minv + 1) {
    i++;
    c = c * 2;
  }
  nBits = i;
}

int DiscreteParameter::getValue() {
  return value;
}

int DiscreteParameter::getNOptions() {
  return maxv - minv + 1;
}

void DiscreteParameter::setBO(int option) {
  value = minv + option;
}
void DiscreteParameter::setHarmonica() {
  // if(std::pow(2,bits.size()) <= maxv-minv+1){
  double v = 0;
  double m = 1;
  for (auto bit : bits) {
    v = v + m * bit.getValue();
    m = m * 2;
  }
  value = static_cast<int>(lround(minv + ((maxv - minv) * (1.0 + v / (m - 1.0)) / 2.0)));
}
} /* namespace datadriven */
} /* namespace sgpp */
