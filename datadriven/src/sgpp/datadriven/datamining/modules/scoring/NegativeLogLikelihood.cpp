// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/scoring/NegativeLogLikelihood.hpp>

#include <iostream>
#include <cmath>

namespace sgpp {
namespace datadriven {

Metric *NegativeLogLikelihood::clone() const { return new NegativeLogLikelihood(*this); }

double NegativeLogLikelihood::measure(const DataVector &predictedValues,
                                      const DataVector &trueValues, const ModelFittingBase &model,
                                      Dataset &testDataset) const {
  DataVector tmp(predictedValues);
  double ll = 0.0;
  for (size_t i = 0; i < predictedValues.size(); i++) {
    double prob = predictedValues.get(i);
    if (prob > 0) {
      ll += std::log(prob);
    }
  }
  return -ll;
}

double NegativeLogLikelihood::measureLowerIsBetter(const DataVector &predictedValues,
                                                   const DataVector &trueValues,
                                                   const ModelFittingBase &model,
                                                   Dataset &testDataset) const {
  return measure(predictedValues, trueValues, model, testDataset);
}
} /* namespace datadriven */
} /* namespace sgpp */
