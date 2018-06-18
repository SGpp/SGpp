/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ContinuousParameter.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: Eric Koepke
 */

#include <cmath>
#include "ContinuousParameter.hpp"

namespace sgpp {
namespace datadriven {

void ContinuousParameter::setHarmonica() {
  double v = 0;
  double m = 1;
  for (auto &bit : bits) {
    v = v + m * bit.getValue();
    m = m * 2;
  }
  value = min + ((max - min) * (1.0 + v / (m - 1.0)) / 2.0);
}

void ContinuousParameter::setBO(double interval) {
  value = (max - min) * interval + min;
}

double ContinuousParameter::getValue() {
  if (logscale) {
    return pow(10, value);
  }
  return value;
}
} /* namespace datadriven */
} /* namespace sgpp */
