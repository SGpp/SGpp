// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSGDE_HPP_
#define LEARNERSGDE_HPP_

#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

struct LearnerSGDEConfiguration {
  // parameters for cross-validation
  bool doCrossValidation_;  // enables cross-validation
  size_t kfold_;  // number of batches for cross validation
  int seed_;  // seed for randomized k-fold
  bool shuffle_;  // randomized/sequential k-fold
  bool silent_;  // verbosity

  // regularization parameter optimization
  double lambda_;  // regularization parameter
  double lambdaStart_;  // lower bound for lambda search range
  double lambdaEnd_;  // upper bound for lambda search range
  // number of lambdas to be tested within the range defined by lambdaStart and lambdaEdns;
  // must be > 1
  size_t lambdaSteps_;
  bool logScale_;  // search the optimization interval on a log-scale
};

// --------------------------------------------------------------------------

class LearnerSGDE: public datadriven::DensityEstimator {
 public:
  /**
   * Constructor
   *
   * @param gridConfig grid configuration
   * @param adaptivityConfig adaptive refinement configuration
   * @param solverConfig solver configuration (CG)
   * @param regularizationConfig config for regularization operator
   * @param learnerSGDEConfig configuration for the learner
   */
  LearnerSGDE(sgpp::base::RegularGridConfiguration& gridConfig,
              sgpp::base::AdpativityConfiguration& adaptivityConfig,
              sgpp::solver::SLESolverConfiguration& solverConfig,
              sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
              LearnerSGDEConfiguration& learnerSGDEConfig);
  virtual ~LearnerSGDE();

  /**
   * Estimate a sparse grid density based on the given data set and
   * the specified configurations.
   *
   * @param samples DataMatrix (nrows = number of samples, ncols = dimensionality)
   */
  void initialize(sgpp::base::DataMatrix& samples);

  /**
   * This methods evaluates the sparse grid density at a single point
   * @param x DataVector length equal to dimensionality
   */
  virtual double pdf(sgpp::base::DataVector& x);

  /**
   * Evaluation of the sparse grid density at a set of points.
   * @param points DataMatrix (nrows = number of samples, ncols = dimensionality)
   * @param res DataVector (size = number of samples) where the results are stored
   */
  virtual void pdf(sgpp::base::DataMatrix& points, sgpp::base::DataVector& res);

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
  virtual void cov(sgpp::base::DataMatrix& cov);

  /**
   * returns the samples in the given dimension
   * @param dim
   */
  virtual sgpp::base::DataVector* getSamples(size_t dim);

  /**
   * returns the complete sample set
   */
  virtual sgpp::base::DataMatrix* getSamples();

  /**
   * get number of dimensions
   */
  virtual size_t getDim();

  /**
   * get number of samples
   */
  virtual size_t getNsamples();

  virtual sgpp::base::Grid* getGrid();
  virtual sgpp::base::DataVector* getAlpha();

  /**
         * returns the surpluses
         */
  virtual sgpp::base::DataVector getSurpluses();

  /**
   * returns the grid storage
   */
  virtual sgpp::base::GridStorage* getGridStorage();

 protected:
  /**
   * Does the learning step on a given grid, training set and regularization parameter lambda
   *
   * @param grid grid
   * @param alpha coefficient vector
   * @param train sample set
   * @param lambdaReg regularization parameter
   */
  virtual void train(sgpp::base::Grid& grid, sgpp::base::DataVector& alpha,
                     sgpp::base::DataMatrix& train, double lambdaReg);

  /**
   * generates a regular grid
   * @param grid grid
   * @param ndim number of dimensions
   */
  void createRegularGrid(sgpp::base::Grid*& grid, size_t ndim);

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
  double computeResidual(sgpp::base::Grid& grid, sgpp::base::DataVector& alpha,
                          sgpp::base::DataMatrix& test, double lambdaReg);

  /**
   * generates the regularization matrix
   * @param grid grid
   */
  sgpp::base::OperationMatrix* computeRegularizationMatrix(
    sgpp::base::Grid& grid);

  /**
   * splits the complete sample set in a set of smaller training and test
   * samples for cross-validation.
   *
   * @param strain vector containing the training samples for cv
   * @param stest vector containing the test samples for cv
   */
  void splitset(std::vector<sgpp::base::DataMatrix*>& strain,
                std::vector<sgpp::base::DataMatrix*>& stest);

  sgpp::base::Grid* grid;
  sgpp::base::DataVector alpha;
  sgpp::base::DataMatrix* samples;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  LearnerSGDEConfiguration learnerSGDEConfig;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERSGDE_HPP_ */

