// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_SCORING_NEGATIVELOGLIKELIHOOD_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_SCORING_NEGATIVELOGLIKELIHOOD_HPP_

#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>

namespace sgpp {
namespace datadriven {
/**
 * Metric that quantifies the likelihood of a dataset given the density function. The smaller
 * the negative log likelihood the better the fit
 */
class NegativeLogLikelihood : public Metric {
 public:
  /**
   * Standard clone method
   * @return the cloned metric instance
   */
  Metric *clone() const override;

  /**
   * Quantify the NLL of the predicted values (i.e. adding the logs of the predicted values
   * and ignorign the true values).
   * Note that probabilities <= 0 are simply ignored (the model can provide those)
   *
   * @param predictedValues probabilites calculated by the model for testing data
   * @param trueValues ignored
   * @param model reference to the model
   * @param testDataset dataset with test data
   * @return the negative log likelihood of the predicted probabilities
   */
  double measure(const DataVector &predictedValues, const DataVector &trueValues,
                 const ModelFittingBase &model, Dataset &testDataset) const override;

  /**
 * Quantify the NLL of the predicted values (i.e. adding the logs of the predicted values
 * and ignorign the true values).
 * Note that probabilities <= 0 are simply ignored (the model can provide those)
 *
 * @param predictedValues probabilites calculated by the model for testing data
 * @param trueValues ignored
 * @param model reference to the model
 * @param testDataset dataset with test data
 * @return the negative log likelihood of the predicted probabilities
 */
  double measureLowerIsBetter(const DataVector &predictedValues, const DataVector &trueValues,
                              const ModelFittingBase &model, Dataset &testDataset) const override;
};
} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_SCORING_NEGATIVELOGLIKELIHOOD_HPP_ */
