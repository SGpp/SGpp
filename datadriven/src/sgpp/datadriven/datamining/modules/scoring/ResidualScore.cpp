// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ResidualScore.hpp>

namespace sgpp {
namespace datadriven {

Metric *ResidualScore::clone() const { return new ResidualScore(*this); }

double ResidualScore::measure(const DataVector &predictedValues, const DataVector &trueValues,
                              const ModelFittingBase &model, Dataset &testDataset) const {
  // get the residual score from the model
  double score = model.computeResidual(testDataset.getData());

  return score;
}

double ResidualScore::measureLowerIsBetter(const DataVector &predictedValues,
                                           const DataVector &trueValues,
                                           const ModelFittingBase &model,
                                           Dataset &testDataset) const {
  return measure(predictedValues, trueValues, model, testDataset);
}

} /* namespace datadriven */
} /* namespace sgpp */
