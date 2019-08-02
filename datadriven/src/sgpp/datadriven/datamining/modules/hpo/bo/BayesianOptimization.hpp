/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * BayesianOptimization.hpp
 *
 *  Created on: Feb 2, 2018
 *      Author: Eric Koepke
 */

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_BAYESIANOPTIMIZATION_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_BAYESIANOPTIMIZATION_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BOConfig.hpp>

#include <vector>
#include "../../../../../../../../base/src/sgpp/base/tools/sle/system/FullSLE.hpp"

namespace sgpp {
namespace datadriven {

/**
 * Class to host all methods to perform Bayesian Optimization
 */
class BayesianOptimization {
 public:
  /**
   * Constructor setting up Gaussian Process
   * @param initialConfigs non-empty vector of initial points to build the Gaussian Process
   */
  explicit BayesianOptimization(const std::vector<BOConfig> &initialConfigs);

  /**
   * Wrapper function for use in optimizer
   * @param inp point in continuous optimization space
   * @return score to optimize on
   */
  double acquisitionOuter(const base::DataVector &inp);

  /**
   * Gaussian Process update step. Incorporates most recent sample into Gaussian Process.
   */
  void updateGP(BOConfig &newConfig, bool normalize);

  /**
   * Implementation of mathematical formulation of the expected improvement acquisition function
   */
  static double acquisitionEI(double dMean, double dVar, double bestsofar);

  /**
   * Perform a Cholesky Decomposition
   * @param km input matrix
   * @param gnew output (triangular) matrix
   */
  void decomposeCholesky(base::DataMatrix &km, base::DataMatrix &gnew);

  /**
   * Solve a system of linear equations using previously decomposed matrix
   * @param gmatrix decomposed matrix
   * @param x target vector
   */
  void solveCholeskySystem(base::DataMatrix &gmatrix, base::DataVector &x);

  /**
   * main routine to find new sample point
   * @param prototype baseline BOConfig
   * @return new sample point
   */
  BOConfig main(BOConfig &prototype);


  /**
   * kernel function
   * @param distance as computed according to some spacial representation
   * @return kernel value representing variance
   */
  double kernel(double distance);

  /**
   * fitting the Gaussian Process to the data in the form of scaling the hyperparameters
   * in relation to each other by likelihood maximization
   */
  base::DataVector fitScales();

  /**
   * Gaussian Process likelihood to perform likelihood maximization over
   * @param inp vector containing scales of the hyperparameter space
   * @return value representative of the likelihood
   */
  double likelihood(const base::DataVector &inp);

  double mean(base::DataVector &knew);

  double var(base::DataVector &knew, double kself);

  void setScales(base::DataVector nscales, double factor);

 protected:
  /**
   * Gram matrix containing all kernel values between all existing samples
   */
  base::DataMatrix kernelmatrix;
  /**
   * Cholesky Decomposition of the Gram matrix
   */
  base::DataMatrix gleft;
  /**
   * solution of Kx = s where K is the Gram matrix and s the scores of the samples
   */
  base::DataVector transformedOutput;
  /**
   * score values belonging to all existing sample points
   */
  base::DataVector rawScores;
  /**
   * scales of the hyperparameter space used to compute the kernel (set by fitScales())
   */
  base::DataVector scales;
  /**
   * best score encountered so far (used for expected improvement calculation)
   */
  double bestsofar;
  /**
   * debugging variable for numerical instabilities
   */
  bool screwedvar;
  /**
   * debugging variable for numerical instabilities
   */
  bool decomFailed = false;
  /**
   * debugging variable for numerical instabilities
   */
  double maxofmax;

  /**
   * existing sample points in the Gaussian Process
   */
  std::vector<BOConfig> allConfigs;
};
} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_BAYESIANOPTIMIZATION_HPP_ */
