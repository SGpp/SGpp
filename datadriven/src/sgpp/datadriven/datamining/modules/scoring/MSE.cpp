// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/scoring/MSE.hpp>

#include <iostream>

namespace sgpp {
namespace datadriven {

Metric *MSE::clone() const { return new MSE(*this); }

double MSE::measure(const DataVector &predictedValues, const DataVector &trueValues,
                    const ModelFittingBase &model, Dataset &testDataset) const {
  DataVector tmp(predictedValues);
  tmp.sub(trueValues);
  // std::cout << "Predicted " << predictedValues.toString() << std::endl;
  // std::cout << "True " << trueValues.toString() << std::endl;
  // std::cout << "Residual " << tmp.toString() << std::endl;
  const double error = tmp.l2Norm();
  return (error * error / static_cast<double>(tmp.getSize()));
}

double MSE::measureLowerIsBetter(const DataVector &predictedValues, const DataVector &trueValues,
                                 const ModelFittingBase &model, Dataset &testDataset) const {
  return measure(predictedValues, trueValues, model, testDataset);
}

} /* namespace datadriven */
} /* namespace sgpp */
