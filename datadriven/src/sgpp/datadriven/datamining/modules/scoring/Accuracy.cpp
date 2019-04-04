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
// FOR EVALUATION ONLY
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/datamining/base/evaluationtools.hpp>

#include <sgpp/datadriven/datamining/modules/scoring/Accuracy.hpp>

#include <iostream>
#include <vector>

using std::vector;

namespace sgpp {
namespace datadriven {

Metric* Accuracy::clone() const { return new Accuracy(*this); }

double Accuracy::measure(const DataVector& predictedValues, const DataVector& trueValues) const {
  size_t correct = 0;

  vector<double> labels{};
  sgpp::base::DataMatrix cM(0, 0, 0.0);

  for (size_t i = 0; i < predictedValues.size(); i++) {
    size_t x = 0;
    size_t y = 0;
    for (size_t a = 0; a < labels.size(); a++) {
      if (trueValues.at(i) == labels.at(a)) {
        x = a + 1;
      }
      if (predictedValues.at(i) == labels.at(a)) {
        y = a + 1;
      }
    }

    // resizing of the classification matrix as response to unseen class labels
    if (x == 0) {
      labels.push_back(trueValues.at(i));
      cM.resizeZero(labels.size(), labels.size());
      x = labels.size();
    }

    if (y == 0) {
      if (trueValues.at(i) != predictedValues.at(i)) {
        labels.push_back(predictedValues.at(i));
        cM.resizeZero(labels.size(), labels.size());
      }
      y = labels.size();
    }
    cM.set(x - 1, y - 1, cM.get(x - 1, y - 1) + 1);
    if (predictedValues.get(i) == trueValues.get(i)) {
      correct++;
    }
  }

  std::cout << evalu::getTime() << "classification matrix:\n" << evalu::dmtS(cM) << std::endl;

  return static_cast<double>(correct) / static_cast<double>(predictedValues.size());
}

} /* namespace datadriven */
} /* namespace sgpp */
