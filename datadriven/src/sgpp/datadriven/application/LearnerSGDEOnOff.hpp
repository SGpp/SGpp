// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#ifndef LEARNERSGDEONOFF_HPP
#define LEARNERSGDEONOFF_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

/**
 * LearnerSGDEOnOff learns the data using sparse grid density estimation. The
 * system matrix is precomputed and factorized using Eigen-, LU- or
 * Cholesky decomposition (offline step). Then, for each class a density
 * function
 * is computed by solving the system in every iteration (online step).
 * If Cholesky decomposition is chosen, refinement/coarsening can be applied.
 */

class LearnerSGDEOnOff {
 public:
  /**
   * Constructor.
   *
   * @param dconf The configuration of the offline object
   * @param trainData The (mandatory) training dataset
   * @param testData The (mandatory) test dataset
   * @param validData The (optional) validation dataset
   * @param classLabels The class labels (e.g. -1, 1)
   * @param classNumber Total number of classes
   * @param usePrior Determines if prior probabilities should be used to compute
   * class labels
   * @param beta The initial weighting factor
   * @param lambda The initial regularization parameter
   */
  LearnerSGDEOnOff(DBMatDensityConfiguration& dconf, Dataset& trainData, Dataset& testData,
                   Dataset* validationData, DataVector& classLabels, size_t classNumber,
                   bool usePrior, double beta, double lambda);

  /**
   * Destructor.
   */
  ~LearnerSGDEOnOff();

  /**
   * Trains the learner with the given dataset.
   *
   * @param batchSize Size of subset of data points used for each training step
   * @param maxDataPasses The number of passes over the whole training data
   * @param refType The refinement indicator (surplus, zero-crossings or
   * data-based)
   * @param refMonitor The refinement strategy (periodic or convergence-based)
   * @param refPeriod The refinement interval (if periodic refinement is chosen)
   * @param accDeclineThreshold The convergence threshold
   *        (if convergence-based refinement is chosen)
   * @param accDeclineBufferSize The number of accuracy measurements which are
   * used to check
   *        convergence (if convergence-based refinement is chosen)
   * @param minRefInterval The minimum number of data points (or data batches)
   * which have to be
   *        processed before next refinement can be scheduled (if
   * convergence-based refinement
   *        is chosen)
   * @param enableCv Specifies whether to perform cross-validation during
   * training process or not
   * @param nextCvStep Determines when next cross-validation has to be triggered
   */
  void train(size_t batchSize, size_t maxDataPasses, std::string refType, std::string refMonitor,
             size_t refPeriod, double accDeclineThreshold, size_t accDeclineBufferSize,
             size_t minRefInterval, bool enableCv, size_t nextCvStep);

  /**
   * Trains the learner with the given data batch
   *
   * @param dataset The next data batch to process
   * @param doCv Enable cross-validation
   * @param refineCoarse Vector of pairs containing a list representing indices
   *        of removed grid points and an unsigned int representing added grid
   * points
   */
  void train(Dataset& dataset, bool doCv = false,
             std::vector<std::pair<std::list<size_t>, size_t> >* refineCoarse = nullptr);

  /**
   * Trains the learner with the given data batch that is already split up wrt
   * its different
   * classes.
   *
   * @param trainDataClasses A vector of pairs; Each pair contains the data
   * points that belong to
   *        one class and the corresponding class label
   * @param doCv Enable cross-validation
   * @param refineCoarse Vector of pairs containing a list representing indices
   * of
   *        removed grid points and an unsigned int representing added grid
   * points
   */
  void train(std::vector<std::pair<DataMatrix*, double> >& trainDataClasses, bool doCv = false,
             std::vector<std::pair<std::list<size_t>, size_t> >* refineCoarse = nullptr);

  /**
   * Returns the accuracy of the classifier measured on the test data.
   *
   * @return The classification accuracy measured on the test data
   */
  double getAccuracy();

  /**
   * Predicts the class labels of the test data points.
   *
   * @param The test data points
   * @return A vector containing the predicted class labels
   */
  DataVector predict(DataMatrix& test);

  /**
   * Predicts the class label of the given data point.
   *
   * @param p The data point
   * @return The predicted class label
   */
  int predict(DataVector& p);

  /**
   * Error evaluation required for convergence-based refinement.
   *
   * @param dataset The data to measure the error on
   * @return The error evaluation
   */
  double getError(Dataset& dataset);

  /**
   * Stores classified data, grids and density function evaluations to csv
   * files.
   */
  void storeResults();

  /**
   * Returns the values of all density functions for a specified data point.
   *
   * @param point The point for which the density functions should be evaluated
   * @return The function evaluations
   */
  DataVector getDensities(DataVector& point);

  /**
   * Sets the cross-validation parameters.
   * They get directly passed to the DBMatOnlineDE class-instance.
   *
   * @param lambdaStep Defines how many different lambdas are tried out
   * @param lambdaStart The smallest possible lambda
   * @param lambdaEnd The biggest possible lambda
   * @param test The test matrix
   * @param testRes The results of the points in the test matrix
   * @param logscale Indicates whether the values between lambdaStart
   *        and lambdaEnd are searched using logscale or not
   */
  void setCrossValidationParameters(int lambdaStep, double lambdaStart, double lambdaEnd,
                                    DataMatrix* test, DataMatrix* testRes, bool logscale);

  /**
  * In case of crossvalidation, returns the current best lambda.
  *
  * @return The lambda value
  */
  // double getBestLambda();

  /**
  * Initialization of online objects in case of Eigen- or LU-decomposition.
  */
  void init();

  /**
  * Returns the number of existing classes.
  *
  * @return The number of classes
  */
  size_t getNumClasses();

  /**
   * Returns the density functions mapped to class labels.
   *
   * @return The density function objects mapped to class labels
   */
  std::vector<std::pair<DBMatOnlineDE*, double> >& getDestFunctions();

  // Stores prior values mapped to class labels
  std::map<double, double> prior;

  // The final classification error
  double error;

  // A vector to store error evaluations
  DataVector avgErrors;

 protected:
  // The training data
  Dataset& trainData;
  // The test data
  Dataset& testData;
  // The (optional) validationData
  Dataset* validationData;

  // The class labels (e.g -1, 1)
  DataVector classLabels;
  // The total number of different classes
  size_t numClasses;

  // Indicates whether the model has been trained or not
  bool trained;
  // Indicates whether the learner has been initialized or not
  bool initDone;
  // Specifies whether prior should be used for class prediction or not
  bool usePrior;
  // Weighting factor
  double beta;

  // The offline object (contains decomposed matrix)
  DBMatOffline offline;
  // The online objects (density functions)
  std::vector<std::pair<DBMatOnlineDE*, double> > destFunctions;

  // Counter for total number of data points processed within ona data pass
  size_t processedPoints;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERSGDEONOFF_HPP */

#endif /* USE_GSL */
