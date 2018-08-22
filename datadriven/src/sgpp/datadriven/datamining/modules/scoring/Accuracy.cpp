/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Accuracy.cpp
 *
 *  Created on: Jul 2, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/modules/scoring/Accuracy.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

Metric* Accuracy::clone() const { return new Accuracy(*this); }

double Accuracy::measure(const DataVector& predictedValues, const DataVector& trueValues) const {
  size_t correct = 0;
  for (size_t i = 0; i < predictedValues.size(); i++) {
    if (predictedValues.get(i) == trueValues.get(i)) {
      correct++;
    }
  }
  return static_cast<double>(correct) / static_cast<double>(predictedValues.size());
}

} /* namespace datadriven */
} /* namespace sgpp */
