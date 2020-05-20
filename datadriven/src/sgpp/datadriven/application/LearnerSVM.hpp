// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/application/PrimalDualSVM.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * LearnerSVM learns the data using support vector machines and sparse grid
 * kernels.
 * As learning algorithm the Pegasos-method is implemented.
 */

class LearnerSVM {
 protected:
  std::unique_ptr<base::Grid> grid;

  base::DataMatrix& trainData;
  base::DataVector& trainLabels;
  base::DataMatrix& testData;
  base::DataVector& testLabels;
  base::DataMatrix* validData;
  base::DataVector* validLabels;

  base::RegularGridConfiguration gridConfig;
  base::AdaptivityConfiguration adaptivityConfig;

  // the svm object
  std::unique_ptr<PrimalDualSVM> svm;

  /**
   * Generates a regular sparse grid.
   *
   * @return The created grid
   */
  std::unique_ptr<base::Grid> createRegularGrid();

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
   * @param pValidData The validation dataset
   * @param pValidLabels The corresponding validation labels
   */
  LearnerSVM(base::RegularGridConfiguration& gridConfig,
             base::AdaptivityConfiguration& adaptivityConfig, base::DataMatrix& pTrainData,
             base::DataVector& pTrainLabels, base::DataMatrix& pTestData,
             base::DataVector& pTestLabels, base::DataMatrix* pValidData,
             base::DataVector* pValidLabels);

  /**
   * Destructor.
   */
  ~LearnerSVM();

  /**
   * Initializes the SVM learner.
   *
   * @param budget The max. number of stored support vectors
   */
  void initialize(size_t budget);

  /**
   * Implements support vector learning with sparse grid kernels.
   *
   * @param maxDataPasses The number of passes over the whole training data
   * @param lambda The regularization parameter
   * @param betaRef Weighting factor for grid points; used within
   *        combined-measure refinement
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

  void train(size_t maxDataPasses, double lambda, double betaRef, std::string refType,
             std::string refMonitor, size_t refPeriod, double errorDeclineThreshold,
             size_t errorDeclineBufferSize, size_t minRefInterval);
  /**
   * Stores classified data, grids and function evaluations to csv files.
   *
   * @param testDataset Data points for which the model is evaluated
   */
  void storeResults(sgpp::base::DataMatrix& testDataset);

  /**
   * Computes the classification accuracy on the given dataset.
   *
   * @param testDataset The data for which class labels should be predicted
   * @param referenceLabels The corresponding actual class labels
   * @param threshold The decision threshold (e.g. for class labels -1, 1 ->
   * threshold = 0)
   * @return The resulting accuracy
   */
  double getAccuracy(sgpp::base::DataMatrix& testDataset,
                     const sgpp::base::DataVector& referenceLabels, const double threshold);

  /**
   * Computes the classification accuracy.
   *
   * @param referenceLabels The actual class labels
   * @param threshold The decision threshold (e.g. for class labels -1, 1 ->
   * threshold = 0)
   * @param predictedLabels The predicted class labels
   * @return The resulting accuracy
   */
  double getAccuracy(const sgpp::base::DataVector& referenceLabels, const double threshold,
                     const sgpp::base::DataVector& predictedLabels);

  /**
   * Predicts class labels based on the trained model.
   *
   * @param testData The data for which class labels should be predicted
   * @param predictedLabels The predicted class labels
   */
  void predict(sgpp::base::DataMatrix& testData, sgpp::base::DataVector& predictedLabels);

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

  // The final classification error
  double error;
  // A vector to store error evaluations
  sgpp::base::DataVector avgErrors;
};

}  // namespace datadriven
}  // namespace sgpp
