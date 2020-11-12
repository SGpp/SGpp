// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * LearnerSGD learns the data using stochastic gradient descent.
 */

class LearnerSGD {
 public:
  /**
   * Constructor.
   *
   * @param gridConfig The grid configuration
   * @param adaptivityConfig The refinement configuration
   * @param pTrainData The training dataset
   * @param pTrainLabels The corresponding training labels
   * @param pTestData The test dataset
   * @param pTestLabels The corresponding test labels
   * @param pValData The validation dataset
   * @param pValLabels The corresponding validation labels
   * @param lambda The regularization parameter
   * @param gamma The learning parameter (i.e. step width)
   * @param batchSize The number of data points which are considered
   *        to compute the error contributions for predictive refinement
   * @param useValidData Specifies if validation data should be used
   *        for all error computations
   */
  LearnerSGD(base::RegularGridConfiguration& gridConfig,
             base::AdaptivityConfiguration& adaptivityConfig, base::DataMatrix& pTrainData,
             base::DataVector& pTrainLabels, base::DataMatrix& pTestData,
             base::DataVector& pTestLabels, base::DataMatrix* pValData,
             base::DataVector* pValLabels, double lambda, double gamma, size_t batchSize,
             bool useValidData);

  /**
   * Destructor.
   */
  ~LearnerSGD();

  /**
   * Initializes the SGD learner (creates grid etc.).
   */
  void initialize();

  /**
   * Implements online learning using stochastic gradient descent.
   *
   * @param maxDataPasses The number of passes over the whole training data
   * @param refType The refinement indicator (surplus, zero-crossings or
   * data-based)
   * @param refMonitor The refinement strategy (periodic or convergence-based)
   * @param refPeriod The refinement interval (if periodic refinement is chosen)
   * @param errorDeclineThreshold The convergence threshold
   *        (if convergence-based refinement is chosen)
   * @param errorDeclineBufferSize The number of error measurements which are
   * used to check
   *        convergence (if convergence-based refinement is chosen)
   * @param minRefInterval The minimum number of data points which have to be
   *        processed before next refinement can be scheduled (if
   * convergence-based refinement
   *        is chosen)
   */
  void train(size_t maxDataPasses, std::string refType, std::string refMonitor, size_t refPeriod,
             double errorDeclineThreshold, size_t errorDeclineBufferSize, size_t minRefInterval);

  /**
   * Computes the classification accuracy on the given dataset.
   *
   * @param testData The data for which class labels should be predicted
   * @param testLabels The corresponding actual class labels
   * @param threshold The decision threshold (e.g. for class labels -1, 1 ->
   * threshold = 0)
   * @return The resulting accuracy
   */
  double getAccuracy(sgpp::base::DataMatrix& testData, sgpp::base::DataVector& testLabels,
                     double threshold);

  /**
   * Stores classified data, grids and function evaluations to csv files.
   *
   * @param testDataset Data points for which the model is evaluated
   */
  void storeResults(base::DataMatrix& testDataset);

  // The final classification error
  double error;
  // A vector to store error evaluations
  sgpp::base::DataVector avgErrors;

 protected:
  /**
   * Generates a regular grid.
   *
   * @return The created grid
   */
  std::unique_ptr<base::Grid> createRegularGrid();

  /**
   * Computes specified error type (e.g. MSE).
   *
   * @param data The data points
   * @param labels The corresponding class labels
   * @param errorType The type of the error measurement (MSE or Hinge loss)
   * @return The error estimation
   */
  double getError(sgpp::base::DataMatrix& data, sgpp::base::DataVector& labels,
                  std::string errorType);

  /**
   * Computes error contribution for each data point of the given
   * data set (required for predictive refinement indicator).
   *
   * @param data The data points
   * @param labels The corresponding class labels
   */
  void getBatchError(sgpp::base::DataMatrix& data, const sgpp::base::DataVector& labels);

  /**
   * Computes the classification accuracy.
   *
   * @param testLabels The actual class labels
   * @param threshold The decision threshold (e.g. for class labels -1, 1 ->
   * threshold = 0)
   * @param predictedLabels The predicted class labels
   * @return The resulting accuracy
   */
  double getAccuracy(sgpp::base::DataVector& testLabels, double threshold,
                     sgpp::base::DataVector& predictedLabels);

  /**
   * Predicts class labels based on the trained model.
   *
   * @param testData The data for which class labels should be predicted
   * @param predictedLabels The predicted class labels
   */
  void predict(base::DataMatrix& testData, base::DataVector& predictedLabels);

  /**
   * Stores the last 'batchSize' processed data points
   * if no validation data is provided.
   *
   * @param x The current data point
   * @param y The corresponding class label
   */
  void pushToBatch(sgpp::base::DataVector& x, double y);

  std::unique_ptr<base::Grid> grid;
  base::DataVector alpha;
  base::DataVector alphaAvg;
  base::DataMatrix& trainData;
  base::DataVector& trainLabels;
  base::DataMatrix& testData;
  base::DataVector& testLabels;
  base::DataMatrix* batchData;
  base::DataVector* batchLabels;
  base::DataVector batchError;

  base::RegularGridConfiguration gridConfig;
  base::AdaptivityConfiguration adaptivityConfig;

  double lambda;
  double gamma;
  double currentGamma;

  size_t batchSize;

  bool useValidData;
};

}  // namespace datadriven
}  // namespace sgpp
