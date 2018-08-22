// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <sgpp/globaldef.hpp>

#include <string>
#include <vector>

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
  // number of lambdas to be tested within the range defined by lambdaStart and lambdaEdns;
  // must be > 1
  size_t lambdaSteps_;
  bool logScale_;  // search the optimization interval on a log-scale
};

struct MakePositiveConfiguration {};

struct SGDEConfiguration {
  bool makePositive_;  // force the density to be positive
  datadriven::MakePositiveCandidateSearchAlgorithm makePositive_candidateSearchAlgorithm_;
  datadriven::MakePositiveInterpolationAlgorithm makePositive_interpolationAlgorithm_;
  bool makePositive_generateConsistentGrid_;
  bool makePositive_verbose_;
  bool unitIntegrand_;  // force unit integrand
};

// --------------------------------------------------------------------------
class SparseGridDensityEstimator;

class SparseGridDensityEstimatorConfiguration : public json::JSON {
  friend class SparseGridDensityEstimator;

 public:
  SparseGridDensityEstimatorConfiguration();
  explicit SparseGridDensityEstimatorConfiguration(const std::string& fileName);

  SparseGridDensityEstimatorConfiguration* clone() override;

  void initConfig();
  sgpp::datadriven::RegularizationType stringToRegularizationType(std::string& regularizationType);
  sgpp::solver::SLESolverType stringToSolverType(std::string& solverType);

 private:
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  sgpp::datadriven::CrossvalidationForRegularizationConfiguration crossvalidationConfig;
  sgpp::datadriven::SGDEConfiguration sgdeConfig;
};

class SparseGridDensityEstimator : public datadriven::DensityEstimator {
 public:
  /**
   * Constructor
   *
   * @param gridConfig grid configuration
   * @param adaptivityConfig adaptive refinement configuration
   * @param solverConfig solver configuration (CG)
   * @param regularizationConfig config for regularization operator
   * @param crossvalidationConfig configuration for the cross validation
   * @param sgdeConfig configuration for the sparse grid density estimation
   */
  SparseGridDensityEstimator(sgpp::base::RegularGridConfiguration& gridConfig,
                             sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                             sgpp::solver::SLESolverConfiguration& solverConfig,
                             sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                             CrossvalidationForRegularizationConfiguration& crossvalidationConfig,
                             SGDEConfiguration& sgdeConfig);

  explicit SparseGridDensityEstimator(SparseGridDensityEstimatorConfiguration& learnerSGDEConfig);

  /**
   * construct learner from given grid and coefficients
   *
   * @param grid
   * @param alpha
   * @param samples
   */
  SparseGridDensityEstimator(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& samples);

  SparseGridDensityEstimator(const SparseGridDensityEstimator& learnerSGDE);

  virtual ~SparseGridDensityEstimator();

  /**
   * Estimate a sparse grid density based on the given data set and
   * the specified configurations.
   *
   * @param samples DataMatrix (nrows = number of samples, ncols = dimensionality)
   */
  virtual void initialize(base::DataMatrix& samples);

  /**
   * This methods evaluates the sparse grid density at a single point
   * @param x DataVector length equal to dimensionality
   */
  double pdf(base::DataVector& x) override;

  /**
   * Evaluation of the sparse grid density at a set of points.
   * @param points DataMatrix (nrows = number of samples, ncols = dimensionality)
   * @param res DataVector (size = number of samples) where the results are stored
   */
  void pdf(base::DataMatrix& points, base::DataVector& res) override;

  /**
   * This method computes the mean of the density function
   */
  double mean() override;

  /**
   * Computes the variance of the density function
   */
  double variance() override;

  /**
   * Computes the covariance of the density function
   */
  void cov(base::DataMatrix& cov, base::DataMatrix* bounds = nullptr) override;

  /**
   * returns the samples in the given dimension
   * @param dim
   */
  std::shared_ptr<base::DataVector> getSamples(size_t dim) override;

  /**
   * returns the complete sample set
   */
  std::shared_ptr<base::DataMatrix> getSamples() override;

  /**
   * get number of dimensions
   */
  size_t getDim() override;

  /**
   * get number of samples
   */
  size_t getNsamples() override;

  /**
   * returns the surpluses
   */
  virtual base::DataVector& getSurpluses();

  /**
   * returns the grid
   */
  virtual base::Grid& getGrid();

  /**
   * Does the learning step on a given grid, training set and regularization parameter lambda
   *
   * @param grid grid
   * @param alpha coefficient vector
   * @param train sample set
   * @param lambdaReg regularization parameter
   */
  virtual void train(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& train,
                     double lambdaReg);

  /**
   * Compute marginal density
   *
   * @param idim dimension which should be left after marginalization
   */
  virtual sgpp::datadriven::SparseGridDensityEstimator* margToDimX(size_t idim);

  /**
   * Marginalize the density in one dimension and return result
   *
   * @param idim dimension which should be marginalized
   */
  virtual sgpp::datadriven::SparseGridDensityEstimator* marginalize(size_t idim);

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
   * @param lambdaReg regularization parameter
   * @return
   */
  double computeResidual(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& test,
                         double lambdaReg);

  /**
   * generates the L^2 dot product matrix
   * @param grid grid
   */
  base::OperationMatrix* computeLTwoDotProductMatrix(base::Grid& grid);

  /**
   * generates the MultipleEval matrix
   * @param grid grid
   * @param train training data
   */
  base::OperationMultipleEval* computeMultipleEvalMatrix(base::Grid& grid, base::DataMatrix& train);

  /**
   * generates the regularization matrix
   * @param grid grid
   */
  base::OperationMatrix* computeRegularizationMatrix(base::Grid& grid);

  /**
   * generates the density system matrix
   * @param grid grid
   * @param train train
   * @param lambdaReg regularization parameter
   */
  std::shared_ptr<datadriven::DensitySystemMatrix> computeDensitySystemMatrix(
      base::Grid& grid, base::DataMatrix& train, double lambdaReg);

  /**
   * splits the complete sample set in a set of smaller training and test
   * samples for cross-validation.
   *
   * @param strain vector containing the training samples for cv
   * @param stest vector containing the test samples for cv
   */
  void splitset(std::vector<std::shared_ptr<base::DataMatrix> >& strain,
                std::vector<std::shared_ptr<base::DataMatrix> >& stest);

  double variance(base::Grid& grid, base::DataVector& alpha);
  double mean(base::Grid& grid, base::DataVector& alpha);

  std::shared_ptr<base::Grid> grid;
  std::shared_ptr<base::DataVector> alpha;

  std::shared_ptr<base::DataMatrix> samples;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  CrossvalidationForRegularizationConfiguration crossvalidationConfig;
  SGDEConfiguration sgdeConfig;
};

}  // namespace datadriven
}  // namespace sgpp
