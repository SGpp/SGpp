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
  vector<vector<size_t>> classificationMatrix{};

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

    if (x == 0) {
      std::cout << "x==0";
      classificationMatrix.resize(classificationMatrix.size() + 1);
      for (size_t g = 0; g < classificationMatrix.size(); g++) {
        std::cout << "x==0 for";
        classificationMatrix.at(g).resize(classificationMatrix.size());
      }
      labels.push_back(trueValues.at(i));
      x = labels.size();
    }

    if (y == 0) {
      std::cout << "y==0";
      if (trueValues.at(i) != predictedValues.at(i)) {
        classificationMatrix.resize(classificationMatrix.size() + 1);
        for (size_t g = 0; g < classificationMatrix.size(); g++) {
          std::cout << "y==0 for";
          classificationMatrix.at(g).resize(classificationMatrix.size());
        }
        labels.push_back(predictedValues.at(i));
      }
      y = labels.size();
    }
    classificationMatrix.at(x - 1).at(y - 1)++;

    std::cout << evalu::getTime() << "classification matrix:\n"
              << evalu::mtS(classificationMatrix) << std::endl;

    if (predictedValues.get(i) == trueValues.get(i)) {
      correct++;
    }
  }
  return static_cast<double>(correct) / static_cast<double>(predictedValues.size());
}

} /* namespace datadriven */
} /* namespace sgpp */
