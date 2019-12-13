// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Metric that quantifies the difference between predicted values and actual values in terms of mean
 * squared error (MSE). MSE is defined strictly positive such that smaller values are better.
 */
class Accuracy : public Metric {
 public:
  Metric* clone() const override;

  /**
   * Quantify the difference between predicted values and actual values in terms of accuracy.
   *
   * @param predictedValues values calculated by the model for testing data
   * @param trueValues actual values as taken from the dataset.
   * @param model reference to the model
   * @param testDataset dataset with test data
   * @return Accuracy - larger values are better.
   */
  double measure(const DataVector& predictedValues, const DataVector& trueValues,
                 const ModelFittingBase& model, Dataset& testDataset) const override;

  /**
   * Quantify the difference between predicted values and actual values in terms of accuracy.
   *
   * @param predictedValues values calculated by the model for testing data
   * @param trueValues actual values as taken from the dataset.
   * @param model reference to the model
   * @param testDataset dataset with test data
   * @return inverse accuracy - smaller values are better.
   */
  double measureLowerIsBetter(const DataVector& predictedValues, const DataVector& trueValues,
                              const ModelFittingBase& model, Dataset& testDataset) const override;
};

} /* namespace datadriven */
} /* namespace sgpp */
