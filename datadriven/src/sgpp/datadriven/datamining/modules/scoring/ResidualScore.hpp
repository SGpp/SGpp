// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Metric that is based on a residual score of the model. For density estimation, this is
 * || R * alpha_lambda - b_val ||_2
 */
class ResidualScore : public Metric {
 public:
  Metric *clone() const override;

  /**
   * Measure the quality of the trained model. Gives the metric access to the trained model, as this
   * is required for some scores.
   *
   * @param predictedValues ignored
   * @param trueValues ignored
   * @param model reference to the model
   * @param testDataset dataset with test data
   */
  double measure(const DataVector &predictedValues, const DataVector &trueValues,
                 const ModelFittingBase &model, Dataset &testDataset) const override;

  /**
   * Measure the quality of the trained model. Gives the metric access to the trained model, as this
   * is required for some scores.
   *
   * @param predictedValues ignored
   * @param trueValues ignored
   * @param model reference to the model
   * @param testDataset dataset with test data
   */
  double measureLowerIsBetter(const DataVector &predictedValues, const DataVector &trueValues,
                              const ModelFittingBase &model, Dataset &testDataset) const override;
};

} /* namespace datadriven */
} /* namespace sgpp */
