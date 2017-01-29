// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSGDE_HPP_
#define LEARNERSGDE_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/datadriven/application/DensityEstimator.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <sgpp/globaldef.hpp>

#include <string>
#include <vector>
#include <map>

namespace sgpp {
namespace datadriven {

struct CrossvalidationForRegularizationConfiguration {
  // parameters for cross-validation
  bool enable_;   // enables cross-validation
  size_t kfold_;  // number of batches for cross validation
  int seed_;      // seed for randomized k-fold
  bool shuffle_;  // randomized/sequential k-fold
  bool silent_;   // verbosity

  // regularization parameter optimization
  double lambda_;       // regularization parameter
  double lambdaStart_;  // lower bound for lambda search range
  double lambdaEnd_;    // upper bound for lambda search range
  // number of lambdas to be tested within the range defined by lambdaStart and
  // lambdaEdns;
  // must be > 1
  size_t lambdaSteps_;
  bool logScale_;  // search the optimization interval on a log-scale
};

// --------------------------------------------------------------------------
class LearnerSGDE;

class LearnerSGDEConfiguration : public json::JSON {
  friend class LearnerSGDE;

 public:
  LearnerSGDEConfiguration();
  explicit LearnerSGDEConfiguration(const std::string& fileName);

  LearnerSGDEConfiguration* clone() override;

  void initConfig();
  sgpp::base::GridType stringToGridType(std::string& gridType);
  sgpp::datadriven::RegularizationType stringToRegularizationType(
      std::string& regularizationType);
  sgpp::solver::SLESolverType stringToSolverType(std::string& solverType);

 private:
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  sgpp::datadriven::CrossvalidationForRegularizationConfiguration
      crossvalidationConfig;
};

class LearnerSGDE : public datadriven::DensityEstimator {
 public:
  /**
   * Constructor
   *
   * @param gridConfig grid configuration
   * @param adaptivityConfig adaptive refinement configuration
   * @param solverConfig solver configuration (CG)
   * @param regularizationConfig config for regularization operator
   * @param crossvalidationConfig configuration for the cross validation
   */
  LearnerSGDE(
      sgpp::base::RegularGridConfiguration& gridConfig,
      sgpp::base::AdpativityConfiguration& adaptivityConfig,
      sgpp::solver::SLESolverConfiguration& solverConfig,
      sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
      CrossvalidationForRegularizationConfiguration& crossvalidationConfig);

  explicit LearnerSGDE(LearnerSGDEConfiguration& learnerSGDEConfig);

  LearnerSGDE(const LearnerSGDE& learnerSGDE);

  virtual ~LearnerSGDE();

  /**
   * Create grid and perform cross-validation if enabled.
   *
   * @param samples DataMatrix (nrows = number of samples, ncols =
   * dimensionality)
   */
  virtual void initialize(base::DataMatrix& samples);

  /**
   * This methods evaluates the sparse grid density at a single point
   * @param x DataVector length equal to dimensionality
   */
  virtual double pdf(base::DataVector& x);

  /**
   * Evaluation of the sparse grid density at a set of points.
   * @param points DataMatrix (nrows = number of samples, ncols =
   * dimensionality)
   * @param res DataVector (size = number of samples) where the results are
   * stored
   */
  virtual void pdf(base::DataMatrix& points, base::DataVector& res);

  /**
   * This method computes the mean of the density function
   */
  virtual double mean();

  /**
   * Computes the variance of the density function
   */
  virtual double variance();

  /**
   * WARNING: Not yet implemented
   */
  virtual void cov(base::DataMatrix& cov);

  /**
   * returns the samples in the given dimension
   * @param dim
   */
  virtual std::shared_ptr<base::DataVector> getSamples(size_t dim);

  /**
   * returns the complete sample set
   */
  virtual std::shared_ptr<base::DataMatrix> getSamples();

  /**
   * get number of dimensions
   */
  virtual size_t getDim();

  /**
   * get number of samples
   */
  virtual size_t getNsamples();

  /**
  * returns the surpluses
  */
  virtual std::shared_ptr<base::DataVector> getSurpluses();

  /**
  * returns the grid
  */
  virtual std::shared_ptr<base::Grid> getGrid();

  /**
   * Does the learning step (i.e. computes pdf) on a given grid,
   * training set and regularization parameter lambda
   *
   * @param grid grid
   * @param alpha coefficient vector
   * @param trainData sample set
   * @param lambdaReg regularization parameter
   */
  virtual void train(base::Grid& grid, base::DataVector& alpha,
                     base::DataMatrix& trainData, double lambdaReg);

  /**
   * Learns the data.
   */
  virtual void train();

  /**
   * Performs the sparse grid density estimation via online learning.
   *
   * @param labels The training labels
   * @param testData The test data
   * @param testLabels The corresponding test labels
   * @param validData The validation data
   * @param validLabels The corresponding validation labels
   * @param classLabels The ocurring class labels (e.g. -1,1)
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
   * @param usePrior Specifies whether prior probabilities should be used to
   * predict class labels
   */
  virtual void trainOnline(base::DataVector& labels, base::DataMatrix& testData,
                           base::DataVector& testLabels, base::DataMatrix* validData,
                           base::DataVector* validLabels, base::DataVector& classLabels,
                           size_t maxDataPasses, std::string refType, std::string refMonitor,
                           size_t refPeriod, double accDeclineThreshold,
                           size_t accDeclineBufferSize, size_t minRefInterval, bool usePrior);

  /**
   * Stores classified data, grids and density function evaluations to csv
   * files.
   *
   * @param testDataset The data for which class labels should be predicted
   */
  virtual void storeResults(base::DataMatrix& testDataset);

  /**
   * Predicts class labels based on the trained model.
   *
   * @param testDataset The data for which class labels should be predicted
   * @param predictedLabels The predicted class labels
   */
  virtual void predict(base::DataMatrix& testDataset,
                       base::DataVector& predictedLabels);

  /**
   * Computes the classification accuracy.
   *
   * @param testDataset The data for which class labels should be predicted
   * @param referenceLabels The corresponding actual class labels
   * @param threshold The decision threshold (e.g. for class labels -1, 1 ->
   * threshold = 0)
   * @return The resulting accuracy
   */
  virtual double getAccuracy(base::DataMatrix& testDataset,
                             const base::DataVector& referenceLabels,
                             const double threshold);

  /**
   * Computes the classification accuracy.
   *
   * @param referenceLabels The actual class labels
   * @param threshold The decision threshold (e.g. for class labels -1, 1 ->
   * threshold = 0)
   * @param predictedLabels The predicted class labels
   * @return The resulting accuracy
   */
  virtual double getAccuracy(const base::DataVector& referenceLabels,
                             const double threshold,
                             const base::DataVector& predictedLabels);

  /**
   * Error evaluation required for convergence-based refinement.
   *
   * @param data The data points to measure the error on
   * @param labels The corresponding class labels
   * @param threshold The decision threshold (e.g. for class labels -1, 1 ->
   * threshold = 0)
   * @param errorType The error type (only "Acc" possible, i.e. classification
   * error
   *        based on accuracy)
   * @return The error evaluation
   */
  virtual double getError(base::DataMatrix& data,
                          const base::DataVector& labels,
                          const double threshold, std::string errorType);

  // The final classification error
  double error;

  // A vector to store error evaluations
  base::DataVector avgErrors;

 protected:
  /**
   * generates a regular grid
   * @return grid
   */
  std::shared_ptr<base::Grid> createRegularGrid();

  /**
   * Does cross-validation to obtain a suitable regularization parameter
   */
  double optimizeLambdaCV();

  /**
   * Compute the residual for a given test data set on a learned grid
   *
   * $|(A - lambda C) alpha - 1/n B|$
   *
   * This is used as quality criterion for the estimated density.
   *
   * @param grid grid
   * @param alpha coefficient vector
   * @param test test set
   * @param lambdaReg regularization parameters
   * @return
   */
  double computeResidual(base::Grid& grid, base::DataVector& alpha,
                         base::DataMatrix& test, double lambdaReg);

  /**
   * generates the regularization matrix
   * @param grid grid
   */
  std::unique_ptr<base::OperationMatrix> computeRegularizationMatrix(
      base::Grid& grid);

  /**
   * splits the complete sample set in a set of smaller training and test
   * samples for cross-validation.
   *
   * @param strain vector containing the training samples for cv
   * @param stest vector containing the test samples for cv
   */
  void splitset(std::vector<std::shared_ptr<base::DataMatrix>>& strain,
                std::vector<std::shared_ptr<base::DataMatrix>>& stest);

  double variance(base::Grid& grid, base::DataVector& alpha);
  double mean(base::Grid& grid, base::DataVector& alpha);

  // the sparse grid
  std::shared_ptr<base::Grid> grid;
  // the density defining coefficient vector (surpluses)
  std::shared_ptr<base::DataVector> alpha;
  // contains sparse grids mapped to class labels (used in trainOnline(...))
  std::map<int, std::shared_ptr<base::Grid>> grids;
  // contains coefficient vectors mapped to class labels (used in
  // trainOnline(...))
  std::map<int, std::shared_ptr<base::DataVector>> alphas;
  // mapping of number of appeared data points to class labels (used in
  // trainOnline(...))
  std::map<int, size_t> appearances;
  // the training data
  std::shared_ptr<base::DataMatrix> trainData;
  // the corresponding class labels
  std::shared_ptr<base::DataVector> trainLabels;
  // stores prior values mapped to class labels
  std::map<int, double> priors;
  // specifies whether prior probabilities should be used to predict class
  // labels
  bool usePrior;
  // regularization parameter
  double lambdaReg;

  // configurations
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  CrossvalidationForRegularizationConfiguration crossvalidationConfig;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERSGDE_HPP_ */
