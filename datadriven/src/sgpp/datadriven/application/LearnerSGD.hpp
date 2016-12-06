// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSGD_HPP
#define LEARNERSGD_HPP

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
   */
  LearnerSGD(sgpp::base::RegularGridConfiguration& gridConfig,
             sgpp::base::AdpativityConfiguration& adaptivityConfig);

  /**
   * Destructor.
   */
  ~LearnerSGD();

  /**
   * Initializes the SGD learner.
   *
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
  void initialize(sgpp::base::DataMatrix& pTrainData,
                  sgpp::base::DataVector& pTrainLabels,
                  sgpp::base::DataMatrix& pTestData,
                  sgpp::base::DataVector& pTestLabels,
                  std::shared_ptr<sgpp::base::DataMatrix> pValData,
                  std::shared_ptr<sgpp::base::DataVector> pValLabels,
                  double lambda, double gamma, size_t batchSize,
                  bool useValidData);

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
  void train(size_t maxDataPasses, std::string refType, std::string refMonitor,
             size_t refPeriod, double errorDeclineThreshold,
             size_t errorDeclineBufferSize, size_t minRefInterval);

  /**
   * Computes the classification accuracy on the given dataset.
   *
   * @param testData The data for which class labels should be predicted
   * @param testLabels The corresponding actual class labels
   * @param threshold The decision threshold (e.g. for class labels -1, 1 ->
   * threshold = 0)
   * @return The resulting accuracy
   */
  double getAccuracy(sgpp::base::DataMatrix& testData,
                     sgpp::base::DataVector& testLabels, double threshold);

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
  std::shared_ptr<base::Grid> createRegularGrid();

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
   * Computes error contribution for each data point of the given data set
   * (required for predictive refinement indicator).
   *
   * @param data The data points
   * @param labels The corresponding class labels
   * @param error The vector which contains the error contributions
   */
  void getBatchError(sgpp::base::DataMatrix& data,
                     sgpp::base::DataVector& labels,
                     sgpp::base::DataVector& error);

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
  void predict(sgpp::base::DataMatrix& testData,
               sgpp::base::DataVector& predictedLabels);

  /**
   * Stores the last 'batchSize' processed data points
   * if no validation data is provided.
   *
   * @param x The current data point
   * @param y The corresponding class label
   */
  void pushToBatch(sgpp::base::DataVector& x, double y);

  std::shared_ptr<base::Grid> grid;
  std::shared_ptr<base::DataVector> alpha;
  std::shared_ptr<base::DataVector> alphaAvg;
  std::shared_ptr<base::DataMatrix> trainData;
  std::shared_ptr<base::DataVector> trainLabels;
  std::shared_ptr<base::DataMatrix> testData;
  std::shared_ptr<base::DataVector> testLabels;
  std::shared_ptr<base::DataMatrix> batchData;
  std::shared_ptr<base::DataVector> batchLabels;
  std::shared_ptr<base::DataVector> batchError;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;

  double lambda;
  double gamma;
  double currentGamma;

  size_t batchSize;

  bool useValidData;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERSGD_HPP */
